#!/usr/bin/env bash
#
# End-to-end validation of the Linux x86 build environments.
#
# Builds MagTense against env-{312,313,314}-linux.yml, with CUDA on and off, and
# asserts that the CUDA toolchain and the newly added packages are actually
# present in the built artefact. Intended to be run on a dedicated Linux box
# after python/.build/regen-envs-linux.sh has updated the env files.
#
# Usage:
#   bash python/.build/test-envs-linux.sh [312|313|314 ...]
#
# Run from the repository root. Set CONDA_DIR if conda is not auto-detected.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

EXPECTED_CUDA="${EXPECTED_CUDA:-12.9}"
MAKE_ARGS=()
[ -n "${CONDA_DIR:-}" ] && MAKE_ARGS+=("CONDA_DIR=${CONDA_DIR}")

CONDA_BIN="${CONDA_BIN:-$(make -s print-CONDA_BIN "${MAKE_ARGS[@]}")}"
ENV_NAME=magtense-env

VERSIONS=("$@")
[ ${#VERSIONS[@]} -eq 0 ] && VERSIONS=(312 313 314)

fail() { echo "FAIL: $*" >&2; exit 1; }
step() { echo; echo "=== $* ==="; }

for V in "${VERSIONS[@]}"; do
    step "python ${V}: rebuilding environment from env-${V}-linux.yml"

    # Makefile:590 only creates the env when it does not already exist, so an
    # edited yml is silently ignored against a stale magtense-env. Removing it
    # first is what makes this test meaningful rather than a false pass.
    make rm-env "${MAKE_ARGS[@]}" || true
    make clean

    step "python ${V}: building with CUDA"
    make python-interface PY_VERSION="$V" USE_CUDA=1 "${MAKE_ARGS[@]}"

    PYTHON="$(make -s print-PYTHON PY_VERSION="$V" "${MAKE_ARGS[@]}")"

    step "python ${V}: asserting CUDA toolchain moved to ${EXPECTED_CUDA}"
    NVCC_VER="$("$CONDA_BIN" run -n "$ENV_NAME" nvcc --version \
        | sed -n 's/.*release \([0-9.]*\),.*/\1/p')"
    echo "nvcc release: ${NVCC_VER}"
    [ "$NVCC_VER" = "$EXPECTED_CUDA" ] \
        || fail "nvcc reports ${NVCC_VER}, expected ${EXPECTED_CUDA}"

    CUDA_PIN="$("$CONDA_BIN" list -n "$ENV_NAME" '^cuda-version$' | awk '/^cuda-version/{print $2}')"
    echo "conda cuda-version: ${CUDA_PIN}"
    [ "$CUDA_PIN" = "$EXPECTED_CUDA" ] \
        || fail "conda cuda-version is ${CUDA_PIN}, expected ${EXPECTED_CUDA}"

    step "python ${V}: checking the extension links CUDA out of the conda prefix"
    SO="$(ls python/src/magtense/lib/magtensesource.cpython-*-x86_64-linux-gnu.so)"
    echo "extension: ${SO}"
    LDD="$(ldd "$SO")"
    echo "$LDD" | grep -E 'libcublas|libcudart|libcusparse' \
        || fail "extension does not link the CUDA libraries"
    if echo "$LDD" | grep -q 'not found'; then
        echo "$LDD" | grep 'not found' >&2
        fail "unresolved shared libraries in ${SO}"
    fi

    step "python ${V}: importing the added packages"
    # pycairo's import name is cairo.
    "$PYTHON" -c "import h5py, scipy, matplotlib, tqdm, importlib_resources, notebook, ipympl, cairo; print('imports OK')" \
        || fail "one of the added packages failed to import"

    step "python ${V}: import magtense and run the test suite"
    "$PYTHON" -c "import magtense; print('magtense', magtense.__file__)"
    make pytest PY_VERSION="$V" "${MAKE_ARGS[@]}"

    step "python ${V}: rebuilding CPU-only"
    # USE_CUDA changes BUILD_FLAGS, and check-flags hard-errors on a mismatch
    # against the previously built libraries, so a clean is required here.
    make clean
    make python-interface PY_VERSION="$V" USE_CUDA=0 "${MAKE_ARGS[@]}"
    "$PYTHON" -c "import magtense; print('cpu build imports OK')"

    echo
    echo "### python ${V}: PASS"
    LAST_V="$V"
done

step "wheel metadata (python ${LAST_V}, cpu variant)"
# dist_pypi.py does not read the extension module from where the build leaves
# it. It expects it staged under python/<variant>_libs/, the way the CI job
# does it (python-package-conda-new.yml:59-60), and it clears
# src/magtense/lib/*.so before copying the staged file in. So stage it here,
# and keep a copy to put back afterwards - otherwise this step silently breaks
# the editable install that the steps above just validated.
BACKUP="$(mktemp -d)"
cp python/src/magtense/lib/magtensesource*.so "$BACKUP/"
mkdir -p python/cpu_libs
cp python/src/magtense/lib/magtensesource*.so python/cpu_libs/

# f2py names the module from the interpreter's EXT_SUFFIX, while dist_pypi.py
# reconstructs the name as cpython-<ver>-x86_64-linux-gnu. A free-threaded or
# otherwise non-standard interpreter makes those disagree, and dist_pypi.py
# reports it only as "the compile step did not produce it", which sends you
# looking at the build instead of at the name.
EXPECTED="magtensesource.cpython-${LAST_V}-x86_64-linux-gnu.so"
if [ ! -f "python/cpu_libs/${EXPECTED}" ]; then
    echo "The staged extension module is not named what dist_pypi.py expects." >&2
    echo "  expected: python/cpu_libs/${EXPECTED}" >&2
    echo "  actually staged:" >&2
    ls -1 python/cpu_libs/ | sed 's/^/    /' >&2
    echo "  interpreter EXT_SUFFIX: $("$PYTHON" -c \
        'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')" >&2
    fail "extension module name does not match dist_pypi.py's expectation"
fi

# Run through conda so that the bare "python -m build" inside dist_pypi.py
# resolves to the environment interpreter rather than whatever is on PATH.
"$CONDA_BIN" run -n "$ENV_NAME" python python/.build/dist_pypi.py \
    --py_version "$LAST_V" --cu_version cpu --platform linux

cp "$BACKUP"/magtensesource*.so python/src/magtense/lib/
rm -rf "$BACKUP" python/cpu_libs

shopt -s nullglob
WHEELS=(python/dist/*.whl)
[ ${#WHEELS[@]} -gt 0 ] || fail "dist_pypi.py produced no wheel"
for whl in "${WHEELS[@]}"; do
    echo "--- $(basename "$whl")"
    META="$(unzip -p "$whl" '*/METADATA')"
    echo "$META" | grep '^Requires-Dist' || true
    echo "$META" | grep -q '^Requires-Dist: ipympl' \
        || fail "ipympl missing from ${whl} metadata"
    echo "$META" | grep -q '^Requires-Dist: pycairo' \
        || fail "pycairo missing from ${whl} metadata"
    # dist_pypi.py:124 strips nvidia-* for the cpu variant.
    if echo "$META" | grep -q '^Requires-Dist: nvidia-'; then
        fail "cpu wheel ${whl} still declares nvidia-* dependencies"
    fi
done

echo
echo "All requested versions passed."
