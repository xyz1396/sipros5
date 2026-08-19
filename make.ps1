# Native Windows build driver.
#
# Install micromamba:
#   Invoke-Expression ((Invoke-WebRequest -Uri https://micro.mamba.pm/install.ps1 -UseBasicParsing).Content)
#
# Environments are installed and maintained by the user; this script only
# validates and uses them. This follows make.sh's environment split:
#   sipros5-release  CPU environment used by build and package
#   sipros5          GPU-capable environment used by buildConda
#
# CPU release environment:
#   micromamba create -n sipros5-release -c conda-forge --strict-channel-priority `
#     python=3.12 cmake ninja "hdf5=2.*=nompi*" `
#     "pytorch-cpu=2.12.1=cpu_mkl*"
#
# NVIDIA GPU environment (PyTorch 2.12.1, CUDA 12.8, Python 3.12):
#   micromamba create -n sipros5 -c conda-forge --strict-channel-priority `
#     python=3.12 cmake ninja "hdf5=2.*=nompi*" `
#     "pytorch-gpu=2.12.1=*cuda128*" "cuda-version=12.8" `
#     "cuda-cudart-dev=12.8" "libcublas-dev=12.8" `
#     "cuda-nvrtc-dev=12.8" "cuda-nvcc=12.8" "cuda-nvtx-dev=12.8"
#
# Visual Studio 2022 with the "Desktop development with C++" workload must be
# installed separately; conda-forge cannot redistribute the MSVC compiler.
#
# Usage:
#   .\make.ps1 build
#   .\make.ps1 buildConda
#   .\make.ps1 package

[CmdletBinding()]
param(
    [ValidateSet('load', 'clean', 'build', 'buildConda', 'debug', 'make', 'package', 'run')]
    [string]$Command = 'build',
    [switch]$EnvironmentActive
)

$ErrorActionPreference = 'Stop'
$RepoDir = $PSScriptRoot
$ReleaseEnvironmentName = if ($env:SIPROS_RELEASE_ENV_NAME) { $env:SIPROS_RELEASE_ENV_NAME } else { 'sipros5-release' }
$CondaEnvironmentName = if ($env:SIPROS_CONDA_ENV_NAME) { $env:SIPROS_CONDA_ENV_NAME } else { 'sipros5' }
$MambaRootPrefix = if ($env:MAMBA_ROOT_PREFIX) { $env:MAMBA_ROOT_PREFIX } else { Join-Path $env:USERPROFILE 'micromamba' }
$MambaExe = if ($env:MAMBA_EXE) {
    $env:MAMBA_EXE
} else {
    $pathMamba = Get-Command micromamba.exe -ErrorAction SilentlyContinue
    if ($pathMamba) { $pathMamba.Source } else { Join-Path $env:LOCALAPPDATA 'micromamba\micromamba.exe' }
}
$BuildRoot = Join-Path $RepoDir 'build\windows'
$ReleaseBuildDir = Join-Path $BuildRoot 'system'
$CondaBuildDir = Join-Path $BuildRoot 'conda'
$DebugBuildDir = Join-Path $BuildRoot 'conda-debug'

function Clear-WindowsOutputs {
    $buildDirectory = Join-Path $RepoDir 'build'
    $binDirectory = Join-Path $RepoDir 'bin'
    if (Test-Path -LiteralPath $buildDirectory) {
        Remove-Item -LiteralPath $buildDirectory -Recurse -Force
    }
    New-Item -ItemType Directory -Path $buildDirectory | Out-Null
    if (Test-Path -LiteralPath $binDirectory) {
        Remove-Item -LiteralPath $binDirectory -Recurse -Force
    }
    New-Item -ItemType Directory -Path $binDirectory | Out-Null

    foreach ($path in @(
        (Join-Path $RepoDir 'tools\aerith.exe'),
        (Join-Path $RepoDir 'tools\sipros.exe'),
        (Join-Path $RepoDir 'tools\siprosMPI.exe'),
        (Join-Path $RepoDir 'siprosRelease.zip')
    )) {
        if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path -Force }
    }
}

function Import-VisualStudioEnvironment {
    if (Get-Command cl.exe -ErrorAction SilentlyContinue) { return }
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
    if (-not (Test-Path -LiteralPath $vswhere)) {
        throw 'Visual Studio 2022 with the Desktop development with C++ workload is required.'
    }
    $installation = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if (-not $installation) { throw 'Install the Visual Studio 2022 Desktop development with C++ workload.' }
    $devShell = Join-Path $installation 'Common7\Tools\VsDevCmd.bat'
    cmd.exe /s /c "`"$devShell`" -no_logo -arch=x64 -host_arch=x64 && set" | ForEach-Object {
        if ($_ -match '^([^=]+)=(.*)$') {
            [Environment]::SetEnvironmentVariable($Matches[1], $Matches[2], 'Process')
        }
    }
    if (-not (Get-Command cl.exe -ErrorAction SilentlyContinue)) {
        throw 'Visual Studio was found, but its x64 C++ compiler could not be activated.'
    }
}

function Test-NvidiaGpu {
    $nvidiaSmi = Get-Command nvidia-smi.exe -ErrorAction SilentlyContinue
    if (-not $nvidiaSmi) { return $false }
    & $nvidiaSmi.Source --query-gpu=name --format=csv,noheader 2>$null | Out-Null
    return $LASTEXITCODE -eq 0
}

function Test-Environment([string]$Name, [bool]$RequireExactTorch, [bool]$RequireGpu = $false) {
    if (-not (Test-Path -LiteralPath $MambaExe)) { return $false }
    $operator = if ($RequireExactTorch) { '==' } else { '>=' }
    $cpuCheck = if ($RequireExactTorch) { ";assert glob.glob(os.path.join(p,'conda-meta','pytorch-2.12.1-cpu_mkl*.json'))" } else { '' }
    $gpuCheck = if ($RequireGpu) { ";assert torch.version.cuda and torch.version.cuda.startswith('12.8')" } else { '' }
    $check = "import glob,os,re,torch;p=os.environ['CONDA_PREFIX'];v=tuple(map(int,re.match(r'\d+(?:\.\d+)+',torch.__version__).group().split('.')));assert v$operator(2,12,1);assert glob.glob(os.path.join(p,'conda-meta','hdf5-2.*-nompi*.json'))$cpuCheck$gpuCheck"
    & $MambaExe run -n $Name python -c $check 2>$null
    return $LASTEXITCODE -eq 0
}

function Require-ReleaseEnvironment {
    if (-not (Test-Path -LiteralPath $MambaExe)) {
        throw "micromamba was not found at $MambaExe. Use the installation command at the top of make.ps1."
    }
    if (-not (Test-Environment $ReleaseEnvironmentName $true)) {
        throw "Missing or stale '$ReleaseEnvironmentName' CPU environment. Install it with the command at the top of make.ps1."
    }
}

function Require-CondaEnvironment {
    if (-not (Test-Path -LiteralPath $MambaExe)) {
        throw "micromamba was not found at $MambaExe. Use the installation command at the top of make.ps1."
    }
    $useGpu = Test-NvidiaGpu
    if (Test-Environment $CondaEnvironmentName $false $useGpu) {
        if ($useGpu) { Write-Host 'buildConda: using the existing CUDA-enabled environment.' }
        else { Write-Host 'buildConda: using the existing CPU environment (no NVIDIA GPU detected).' }
        return
    }
    if ($useGpu) {
        throw "Missing or stale '$CondaEnvironmentName' GPU environment. Install the CUDA 12.8 environment with the command at the top of make.ps1."
    } else {
        throw "Missing or stale '$CondaEnvironmentName' environment. Install it before running buildConda."
    }
}

function Enter-MicromambaEnvironment([string]$Name) {
    $powershell = (Get-Process -Id $PID).Path
    & $MambaExe run -n $Name $powershell -NoProfile -ExecutionPolicy Bypass `
        -File $PSCommandPath -Command $Command -EnvironmentActive
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}

function Invoke-WindowsBuild([string]$Directory, [string]$Configuration) {
    $torchPrefix = python -c 'import torch; print(torch.utils.cmake_prefix_path)'
    if ($LASTEXITCODE -ne 0) { throw 'Unable to locate the PyTorch CMake package.' }
    $torchDir = Join-Path $torchPrefix.Trim() 'Torch'
    $hdf5Dir = Join-Path $env:CONDA_PREFIX 'Library\cmake'
    $arguments = @(
        '-S', $RepoDir, '-B', $Directory, '-G', 'Ninja',
        "-DCMAKE_BUILD_TYPE=$Configuration", '-DBUILD_CONDA=ON',
        '-DSIPROS_BUILD_MPI=OFF', '-DAERITH_ENABLE_TORCH=ON',
        "-DTorch_DIR=$torchDir", '-DAERITH_TORCH_PACKAGE_VERSION=2.12.1',
        "-DHDF5_DIR=$hdf5Dir", "-DHDF5_ROOT=$env:CONDA_PREFIX",
        '-DHDF5_USE_STATIC_LIBRARIES=OFF',
        '-DCMAKE_DISABLE_FIND_PACKAGE_ZLIB=TRUE',
        '-DOpenMP_C_FLAGS=/openmp:llvm',
        '-DOpenMP_CXX_FLAGS=/openmp:llvm'
    )
    # Only buildConda/debug may opt into CUDA. The CPU release/package build
    # must not inherit an unrelated system CUDA installation from PATH.
    $nvcc = if ($Command -in 'buildConda', 'debug') {
        Get-Command nvcc.exe -ErrorAction SilentlyContinue
    } else {
        $null
    }
    if ($null -ne $nvcc) {
        $arguments += '-DAERITH_ENABLE_CUDA=ON'
        $arguments += "-DCMAKE_CUDA_COMPILER=$($nvcc.Source)"
        $arguments += "-DCUDAToolkit_ROOT=$env:CONDA_PREFIX"
        Write-Host "Configuring CUDA with $($nvcc.Source)"
    } else {
        $arguments += '-DAERITH_ENABLE_CUDA=OFF'
    }
    & cmake @arguments
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    & cmake --build $Directory --target sipros aerith
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}

function Publish-WindowsTools([string]$SourceBinDirectory) {
    $executables = @('sipros.exe', 'aerith.exe')
    $models = @(
        'diann-2.6.1-fragmentation.pt',
        'diann-2.6.1-retention-time.pt'
    )
    foreach ($name in @($executables + $models)) {
        $source = Join-Path $SourceBinDirectory $name
        if (-not (Test-Path -LiteralPath $source)) {
            throw "Missing built runtime file: $source"
        }
    }
    $binDestination = Join-Path $RepoDir 'bin'
    $toolsDestination = Join-Path $RepoDir 'tools'
    New-Item -ItemType Directory -Path $binDestination -Force | Out-Null
    New-Item -ItemType Directory -Path $toolsDestination -Force | Out-Null
    foreach ($name in $models) {
        $staleModel = Join-Path $binDestination $name
        if (Test-Path -LiteralPath $staleModel) {
            Remove-Item -LiteralPath $staleModel -Force
        }
    }
    foreach ($name in $executables) {
        Copy-Item -LiteralPath (Join-Path $SourceBinDirectory $name) `
            -Destination (Join-Path $binDestination $name) -Force
    }
    foreach ($name in @($executables + $models)) {
        Copy-Item -LiteralPath (Join-Path $SourceBinDirectory $name) `
            -Destination (Join-Path $toolsDestination $name) -Force
    }
}

function New-WindowsPackage([string]$SourceBinDirectory) {
    $stageRoot = Join-Path $env:TEMP ("sipros5-package-" + [guid]::NewGuid().ToString('N'))
    $stage = Join-Path $stageRoot 'sipros'
    $archive = Join-Path $RepoDir 'siprosRelease.zip'
    try {
        $tools = Join-Path $stage 'tools'
        New-Item -ItemType Directory -Path $tools -Force | Out-Null
        Copy-Item -LiteralPath (Join-Path $RepoDir 'script33') -Destination $stage -Recurse
        Copy-Item -LiteralPath (Join-Path $RepoDir 'LICENSE') -Destination $stage
        foreach ($name in 'sipros.exe', 'aerith.exe', 'diann-2.6.1-fragmentation.pt', 'diann-2.6.1-retention-time.pt') {
            $source = Join-Path $SourceBinDirectory $name
            if (-not (Test-Path -LiteralPath $source)) { throw "Missing built runtime file: $source" }
            Copy-Item -LiteralPath $source -Destination $tools
        }
        # Windows searches beside the executable for dependent DLLs.
        Get-ChildItem -LiteralPath (Join-Path $env:CONDA_PREFIX 'Library\bin') -Filter '*.dll' | Copy-Item -Destination $tools
        if (Test-Path -LiteralPath $archive) { Remove-Item -LiteralPath $archive -Force }
        Compress-Archive -LiteralPath $stage -DestinationPath $archive -CompressionLevel Optimal
        Write-Host "Package created: $archive"
    } finally {
        if (Test-Path -LiteralPath $stageRoot) { Remove-Item -LiteralPath $stageRoot -Recurse -Force }
    }
}

$env:MAMBA_ROOT_PREFIX = $MambaRootPrefix
if ($Command -eq 'clean') {
    Clear-WindowsOutputs
    return
}
if ($Command -eq 'load') { return }

Import-VisualStudioEnvironment
if (-not $EnvironmentActive) {
    if ($Command -in 'build', 'package', 'make', 'run') {
        Require-ReleaseEnvironment
        Enter-MicromambaEnvironment $ReleaseEnvironmentName
    } elseif ($Command -in 'buildConda', 'debug') {
        Require-CondaEnvironment
        Enter-MicromambaEnvironment $CondaEnvironmentName
    }
    return
}

switch ($Command) {
    'build' {
        Invoke-WindowsBuild $ReleaseBuildDir 'Release'
        Publish-WindowsTools (Join-Path $ReleaseBuildDir 'bin')
        Write-Host "Windows release binaries published in $RepoDir\bin"
    }
    'buildConda' {
        Invoke-WindowsBuild $CondaBuildDir 'Release'
        Publish-WindowsTools (Join-Path $CondaBuildDir 'bin')
        Write-Host "Windows Conda binaries published in $RepoDir\bin"
    }
    'debug' { Invoke-WindowsBuild $DebugBuildDir 'Debug' }
    'make' {
        & cmake --build $ReleaseBuildDir
        if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
    }
    'package' {
        Clear-WindowsOutputs
        Invoke-WindowsBuild $ReleaseBuildDir 'Release'
        Publish-WindowsTools (Join-Path $ReleaseBuildDir 'bin')
        New-WindowsPackage (Join-Path $ReleaseBuildDir 'bin')
    }
    'run' {
        Write-Host 'Use the HDF5 workflow entrypoint, for example:'
        Write-Host '  micromamba run -n sipros5 python script33/main.py <arguments>'
    }
}
