<#
.SYNOPSIS
    One-shot installer for a MagTense development environment on Windows.

.DESCRIPTION
    Provisions everything needed to build and use the MagTense Fortran core and
    Python interface, WITHOUT requiring a manual Visual Studio / Intel oneAPI
    system install. The entire toolchain (Intel ifx/icx, MKL, the MSVC build
    tools, make, ninja, python, numpy) is provided by a pinned conda environment
    ("magtense-env"); only CVODE/SUNDIALS is compiled from source.

    The script is idempotent: re-running it resumes from wherever it left off
    (existing conda install, existing env, already-built CVODE are all detected
    and skipped).

    Steps:
      1. Preflight + start a transcript log.
      2. Obtain the MagTense source (use the surrounding clone, or download a zip).
      3. Install Miniforge3 (userspace) if conda is not already available.
      4. Create the "magtense-env" conda environment (GPU or CPU variant).
      5. Download + cmake-build CVODE 7.4.0 into <repo>\cvode.
      6. Build the Fortran core + Python interface and pip-install it (editable).
      7. Verify "import magtense".
      8. Create a Start Menu shortcut that opens a shell with the env activated.

.PARAMETER InstallDir
    Directory that will hold the MagTense source when it must be downloaded.
    Default: %USERPROFILE%\MagTense. Ignored when the script is run from inside
    an existing clone.

.PARAMETER Compute
    'gpu', 'cpu', or 'auto' (default). 'auto' picks GPU when an NVIDIA adapter is
    detected, otherwise CPU. Selects the conda env file and the make USE_CUDA flag.

.PARAMETER PyVersion
    Python minor version for the env file: '314' (default), '313', or '312'.

.PARAMETER Ref
    Git ref (branch, tag, or SHA) to download when not running from a clone.
    Default: 'master'.

.PARAMETER CondaDir
    Where to install Miniforge3 if conda is not found.
    Default: %USERPROFILE%\miniforge3.

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File install-magtense.ps1

.EXAMPLE
    powershell -ExecutionPolicy Bypass -File install-magtense.ps1 -Compute cpu -PyVersion 313
#>

[CmdletBinding()]
param(
    [string]$InstallDir = (Join-Path $env:USERPROFILE 'MagTense'),
    [ValidateSet('auto', 'gpu', 'cpu')]
    [string]$Compute = 'auto',
    [ValidateSet('314', '313', '312')]
    [string]$PyVersion = '314',
    [string]$Ref = 'master',
    [string]$CondaDir = (Join-Path $env:USERPROFILE 'miniforge3')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

# --- constants -------------------------------------------------------------
$EnvName       = 'magtense-env'
$GitHubOwner   = 'cmt-dtu-energy'
$GitHubRepo    = 'MagTense'
$CvodeVersion  = '7.4.0'
$MiniforgeUrl  = 'https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Windows-x86_64.exe'

# --- helpers ---------------------------------------------------------------
function Write-Step   { param([string]$Msg) Write-Host "`n==> $Msg" -ForegroundColor Cyan }
function Write-Info   { param([string]$Msg) Write-Host "    $Msg" -ForegroundColor Gray }
function Write-Ok     { param([string]$Msg) Write-Host "    $Msg" -ForegroundColor Green }

function Invoke-Native {
    # Run an external command and throw with context if it fails.
    param([Parameter(Mandatory)][scriptblock]$Command, [string]$What = 'command')
    & $Command
    if ($LASTEXITCODE -ne 0) {
        throw "$What failed with exit code $LASTEXITCODE."
    }
}

# --- 1. preflight ----------------------------------------------------------
$LogDir = Join-Path $env:USERPROFILE '.magtense'
New-Item -ItemType Directory -Force -Path $LogDir | Out-Null
$LogFile = Join-Path $LogDir ("install-{0:yyyyMMdd-HHmmss}.log" -f (Get-Date))
Start-Transcript -Path $LogFile -Force | Out-Null

try {
    Write-Step "MagTense Windows installer"
    Write-Info "Log file: $LogFile"

    if (-not [Environment]::Is64BitOperatingSystem) {
        throw "MagTense requires 64-bit Windows."
    }

    # Resolve compute mode
    if ($Compute -eq 'auto') {
        $nvidia = @(Get-CimInstance Win32_VideoController -ErrorAction SilentlyContinue |
                    Where-Object { $_.Name -match 'NVIDIA' })
        if ($nvidia.Count -gt 0) {
            $Compute = 'gpu'
            Write-Info "NVIDIA GPU detected -> building with CUDA (Compute=gpu)."
        } else {
            $Compute = 'cpu'
            Write-Info "No NVIDIA GPU detected -> CPU-only build (Compute=cpu)."
        }
    }
    $UseCuda = if ($Compute -eq 'gpu') { 1 } else { 0 }

    # --- 2. obtain source --------------------------------------------------
    Write-Step "Locating MagTense source"
    $repoFromScript = $null
    if ($PSScriptRoot) {
        $candidate = (Resolve-Path (Join-Path $PSScriptRoot '..\..') -ErrorAction SilentlyContinue)
        if ($candidate -and
            (Test-Path (Join-Path $candidate 'Makefile')) -and
            (Test-Path (Join-Path $candidate 'python\.build'))) {
            $repoFromScript = $candidate.Path
        }
    }

    if ($repoFromScript) {
        $RepoDir = $repoFromScript
        Write-Ok "Using existing checkout: $RepoDir"
    } else {
        New-Item -ItemType Directory -Force -Path $InstallDir | Out-Null
        $RepoDir = Join-Path $InstallDir 'MagTense'
        if (Test-Path (Join-Path $RepoDir 'Makefile')) {
            Write-Ok "Source already present: $RepoDir"
        } else {
            $zipUrl = "https://github.com/$GitHubOwner/$GitHubRepo/archive/$Ref.zip"
            $zipPath = Join-Path $env:TEMP "magtense-$Ref.zip"
            Write-Info "Downloading $zipUrl"
            Invoke-WebRequest -Uri $zipUrl -OutFile $zipPath -UseBasicParsing
            Write-Info "Extracting to $InstallDir"
            Expand-Archive -Path $zipPath -DestinationPath $InstallDir -Force
            # GitHub zips extract to <repo>-<ref-with-slashes-as-dashes>
            $extracted = Get-ChildItem -Path $InstallDir -Directory |
                         Where-Object { $_.Name -like "$GitHubRepo-*" } |
                         Sort-Object LastWriteTime -Descending | Select-Object -First 1
            if (-not $extracted) { throw "Could not find extracted source under $InstallDir." }
            Rename-Item -Path $extracted.FullName -NewName 'MagTense'
            Remove-Item $zipPath -Force
            Write-Ok "Source ready: $RepoDir"
        }
    }

    # --- 3. install Miniforge ---------------------------------------------
    Write-Step "Ensuring conda (Miniforge3) is installed"
    $CondaExe = Join-Path $CondaDir 'Scripts\conda.exe'
    if (-not (Test-Path $CondaExe)) {
        # Maybe conda is already installed from another distribution. Only an
        # Application has a real path in .Source; if `conda init` has run, the
        # PATH entry is a PowerShell *function* whose .Source is empty, so we
        # restrict the lookup to Application and fall back to $env:CONDA_EXE
        # (which conda sets to the real conda.exe path when it is initialized).
        $existing = Get-Command conda -CommandType Application -ErrorAction SilentlyContinue |
                    Select-Object -First 1
        $existingExe = if ($existing) { $existing.Source } elseif ($env:CONDA_EXE) { $env:CONDA_EXE } else { $null }
        if ($existingExe -and (Test-Path $existingExe)) {
            $CondaExe = $existingExe
            # Ask conda for its base dir; the executable can live at various
            # depths (…\Scripts\conda.exe, …\condabin\conda.bat,
            # …\Library\bin\conda.bat) so string-splitting the path is unreliable.
            $base = (& $CondaExe info --base 2>$null | Select-Object -First 1)
            if ($base) { $base = $base.Trim() }
            if ($base -and (Test-Path $base)) {
                $CondaDir = $base
            } else {
                $CondaDir = Split-Path (Split-Path $CondaExe -Parent) -Parent
            }
            Write-Ok "Using existing conda: $CondaExe (base: $CondaDir)"
        } else {
            $installer = Join-Path $env:TEMP 'Miniforge3-Windows-x86_64.exe'
            Write-Info "Downloading Miniforge3..."
            Invoke-WebRequest -Uri $MiniforgeUrl -OutFile $installer -UseBasicParsing
            Write-Info "Installing Miniforge3 to $CondaDir (silent, userspace)..."
            $mfArgs = @('/InstallationType=JustMe', '/RegisterPython=0', '/AddToPath=0', '/S', "/D=$CondaDir")
            $p = Start-Process -FilePath $installer -ArgumentList $mfArgs -Wait -PassThru
            if ($p.ExitCode -ne 0) { throw "Miniforge installer exited with code $($p.ExitCode)." }
            Remove-Item $installer -Force
            $CondaExe = Join-Path $CondaDir 'Scripts\conda.exe'
            if (-not (Test-Path $CondaExe)) { throw "conda.exe not found after install at $CondaExe." }
            Write-Ok "Miniforge3 installed at $CondaDir"
        }
    } else {
        Write-Ok "conda already installed at $CondaDir"
    }

    # Bring conda into this session so we can `conda activate` the build env.
    # Activation is essential: the vs2022/Intel conda packages ship activate.d
    # scripts that set INCLUDE/LIB/PATH for the MSVC + ifx toolchain. Without
    # activation the compilers/linker would not find their environment.
    $CondaHook = Join-Path $CondaDir 'shell\condabin\conda-hook.ps1'
    if (-not (Test-Path $CondaHook)) { throw "conda-hook.ps1 not found at $CondaHook." }
    & $CondaHook

    # --- 4. create env -----------------------------------------------------
    Write-Step "Creating conda environment '$EnvName' ($Compute build)"
    # Resolve the env by *prefix*, not by a substring match on `conda env list`.
    # That text listing includes stale entries from ~/.conda/environments.txt
    # (shared across every conda/miniconda/miniforge on the machine and left
    # behind when an env is deleted or built by a different distribution). A
    # substring match on such a phantom path makes us skip creation, and then
    # `conda activate <name>` fails with EnvironmentNameNotFound because the
    # prefix is not under *this* conda's envs_dirs. Instead, list prefixes as
    # JSON and keep only one whose leaf is $EnvName and that actually contains a
    # conda-meta dir (i.e. a real, populated env), then activate it by path.
    $envJson  = & $CondaExe env list --json | ConvertFrom-Json
    $EnvPrefix = $envJson.envs |
        Where-Object { (Split-Path $_ -Leaf) -eq $EnvName -and
                       (Test-Path (Join-Path $_ 'conda-meta')) } |
        Select-Object -First 1
    if ($EnvPrefix) {
        Write-Ok "Environment '$EnvName' already exists at $EnvPrefix (skipping creation)."
    } else {
        $suffix  = if ($Compute -eq 'cpu') { '-cpu' } else { '' }
        $envFile = Join-Path $RepoDir "python\.build\env-$PyVersion-win$suffix.yml"
        if (-not (Test-Path $envFile)) { throw "Env file not found: $envFile" }
        Write-Info "conda env create -f $envFile"
        Invoke-Native { & $CondaExe env create -n $EnvName -f $envFile } 'conda env create'
        # Re-resolve the prefix conda just created so we activate by path below.
        $envJson  = & $CondaExe env list --json | ConvertFrom-Json
        $EnvPrefix = $envJson.envs |
            Where-Object { (Split-Path $_ -Leaf) -eq $EnvName -and
                           (Test-Path (Join-Path $_ 'conda-meta')) } |
            Select-Object -First 1
        if (-not $EnvPrefix) { throw "Environment '$EnvName' not found after creation." }
        Write-Ok "Environment created at $EnvPrefix."
    }

    # Activate by prefix: name-based activation only works when the prefix's
    # parent is a registered envs_dir, whereas a prefix path always resolves.
    conda activate $EnvPrefix
    Write-Ok "Activated '$EnvName' ($EnvPrefix)."

    # --- 5. build CVODE ----------------------------------------------------
    Write-Step "Building CVODE $CvodeVersion"
    $CvodeRoot = Join-Path $RepoDir 'cvode'
    if (Test-Path (Join-Path $CvodeRoot 'lib')) {
        Write-Ok "CVODE already built at $CvodeRoot (skipping)."
    } else {
        # cmake is required to build CVODE but is not part of the pinned env
        # spec; install it into the env on demand (idempotent).
        if (-not (Get-Command cmake -ErrorAction SilentlyContinue)) {
            Write-Info "cmake not found in env; installing via conda..."
            Invoke-Native { & $CondaExe install -p $EnvPrefix -y cmake } 'conda install cmake'
        }
        $tgz = Join-Path $env:TEMP "cvode-$CvodeVersion.tar.gz"
        $srcParent = Join-Path $env:TEMP "cvode-src-$CvodeVersion"
        $cvUrl = "https://github.com/LLNL/sundials/releases/download/v$CvodeVersion/cvode-$CvodeVersion.tar.gz"
        Write-Info "Downloading $cvUrl"
        Invoke-WebRequest -Uri $cvUrl -OutFile $tgz -UseBasicParsing
        if (Test-Path $srcParent) { Remove-Item $srcParent -Recurse -Force }
        New-Item -ItemType Directory -Force -Path $srcParent | Out-Null
        Write-Info "Extracting CVODE sources"
        # Prefer the native Windows bsdtar (System32\tar.exe): it handles
        # "C:\..." drive-letter paths directly. The env's GNU tar (from the
        # git/MSYS package) instead reads the archive path as a remote
        # host:path spec and mangles the backslash paths it is given.
        $winTar = Join-Path $env:SystemRoot 'System32\tar.exe'
        if (Test-Path $winTar) {
            Invoke-Native { & $winTar -xf $tgz -C $srcParent } 'tar (extract cvode)'
        } else {
            # Fallback for older Windows without bsdtar: GNU tar needs
            # forward-slash paths plus --force-local for the drive colon.
            $tgzFwd = $tgz -replace '\\', '/'
            $srcFwd = $srcParent -replace '\\', '/'
            Invoke-Native { & tar --force-local -xf $tgzFwd -C $srcFwd } 'tar (extract cvode)'
        }
        $CvodeSrc = (Get-ChildItem -Path $srcParent -Directory | Select-Object -First 1).FullName

        $CvodeBuild = Join-Path $env:TEMP "cvode-build-$CvodeVersion"
        if (Test-Path $CvodeBuild) { Remove-Item $CvodeBuild -Recurse -Force }

        # Configure/build/install with the env's Intel compilers + Ninja.
        # Flags mirror python/README.md and .github/workflows/cmake-sundials-cvode.yml.
        Write-Info "Configuring CVODE (cmake -G Ninja, ifx/icx-cl)"
        Invoke-Native {
            & cmake -G Ninja -B $CvodeBuild -S $CvodeSrc `
                -D CMAKE_BUILD_TYPE=Release `
                -D BUILD_ARKODE=OFF -D BUILD_CVODE=ON -D BUILD_CVODES=OFF `
                -D BUILD_IDA=OFF -D BUILD_IDAS=OFF -D BUILD_KINSOL=OFF `
                -D BUILD_SHARED_LIBS=OFF -D BUILD_STATIC_LIBS=ON `
                -D CMAKE_INSTALL_PREFIX=$CvodeRoot `
                -D EXAMPLES_INSTALL_PATH="$CvodeRoot\examples" `
                -D CMAKE_C_FLAGS=-Wno-deprecated-declarations `
                -D CMAKE_C_COMPILER=icx-cl -D CMAKE_CXX_COMPILER=icx-cl `
                -D CMAKE_Fortran_COMPILER=ifx `
                -D BUILD_FORTRAN_MODULE_INTERFACE=ON -D ENABLE_OPENMP=ON
        } 'cmake configure (cvode)'
        Write-Info "Building + installing CVODE"
        Invoke-Native { & cmake --build $CvodeBuild --config Release } 'cmake --build (cvode)'
        Invoke-Native { & cmake --install $CvodeBuild } 'cmake --install (cvode)'
        Remove-Item $tgz -Force
        Write-Ok "CVODE installed to $CvodeRoot"
    }

    # --- 6. build core + interface ----------------------------------------
    Write-Step "Building Fortran core + Python interface (USE_CUDA=$UseCuda)"
    Push-Location $RepoDir
    try {
        Invoke-Native { & make python-interface-win USE_CUDA=$UseCuda } 'make python-interface-win'
    } finally {
        Pop-Location
    }
    Write-Ok "Core + interface built and installed (editable) into '$EnvName'."

    # --- 7. verify ---------------------------------------------------------
    Write-Step "Verifying installation"
    Invoke-Native { & python -c "import magtense; print('magtense OK:', magtense.__file__)" } 'import magtense'
    Write-Ok "import magtense succeeded."

    # --- 8. Start Menu shortcut -------------------------------------------
    Write-Step "Creating Start Menu shortcut"
    try {
        $programs = [Environment]::GetFolderPath('Programs')
        $lnk = Join-Path $programs 'MagTense Dev Shell.lnk'
        $ws = New-Object -ComObject WScript.Shell
        $sc = $ws.CreateShortcut($lnk)
        $sc.TargetPath = "$env:SystemRoot\System32\WindowsPowerShell\v1.0\powershell.exe"
        $sc.Arguments  = "-NoExit -ExecutionPolicy Bypass -Command " +
                         "`"& '$CondaHook'; conda activate $EnvName; Set-Location '$RepoDir'`""
        $sc.WorkingDirectory = $RepoDir
        $sc.Description = 'PowerShell with the MagTense conda environment activated'
        $sc.Save()
        Write-Ok "Shortcut created: $lnk"
    } catch {
        Write-Info "Could not create Start Menu shortcut ($($_.Exception.Message)). Skipping."
    }

    Write-Step "Done!"
    Write-Host ""
    Write-Host "MagTense is installed. To start working:" -ForegroundColor Green
    Write-Host "  * Open 'MagTense Dev Shell' from the Start Menu, or run:" -ForegroundColor Green
    Write-Host "      & '$CondaHook'; conda activate $EnvName" -ForegroundColor Green
    Write-Host "  * Repository:      $RepoDir" -ForegroundColor Green
    Write-Host "  * Example scripts: $RepoDir\python\examples" -ForegroundColor Green
    Write-Host ""
}
catch {
    Write-Host ""
    Write-Host "INSTALL FAILED: $($_.Exception.Message)" -ForegroundColor Red
    Write-Host "See the full log at: $LogFile" -ForegroundColor Red
    Write-Host "You can safely re-run this installer; completed steps are skipped." -ForegroundColor Yellow
    exit 1
}
finally {
    Stop-Transcript | Out-Null
}
