#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
TEST_ROOT="${SCRIPT_DIR}"

if ! command -v fpm >/dev/null 2>&1; then
    echo "ERROR: fpm not found on PATH."
    echo "Install with: conda install -c conda-forge fpm"
    exit 1
fi

found_any=0

for manifest in "${TEST_ROOT}"/*/fpm.toml; do
    if [[ -f "${manifest}" ]]; then
        found_any=1
        harness_dir="$(cd -- "$(dirname -- "${manifest}")" && pwd)"
        harness_name="$(basename -- "${harness_dir}")"

        echo "================================================================================"
        echo "Running fpm tests: ${harness_name}"
        echo "Directory: ${harness_dir}"
        echo "================================================================================"

        ( cd "${harness_dir}" && fpm test )
    fi
done

if [[ "${found_any}" -eq 0 ]]; then
    echo "ERROR: No test harnesses found under '${TEST_ROOT}/*/fpm.toml'."
    exit 1
fi
