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
#     "pytorch-cpu=2.12.1=cpu_mkl*" imgui=1.92.9 libvulkan-headers `
#     "dotnet-sdk=8.*"
#
# NVIDIA GPU environment (PyTorch 2.12.1, CUDA 12.8, Python 3.12):
#   micromamba create -n sipros5 -c conda-forge --strict-channel-priority `
#     python=3.12 cmake ninja "hdf5=2.*=nompi*" `
#     "pytorch-gpu=2.12.1=*cuda128*" "cuda-version=12.8" `
#     "cuda-cudart-dev=12.8" "libcublas-dev=12.8" `
#     "cuda-nvrtc-dev=12.8" "cuda-nvcc=12.8" "cuda-nvtx-dev=12.8" `
#     imgui=1.92.9 libvulkan-headers
#
# Visual Studio 2022 with the "Desktop development with C++" workload must be
# installed separately; conda-forge cannot redistribute the MSVC compiler.
#
# Usage:
#   .\make.ps1 build
#   .\make.ps1 buildConda
#   .\make.ps1 package  # creates sipros_windows_<version>.zip and .msi
# Set SIPROS_VERSION to override the default three-part MSI version.

[CmdletBinding()]
param(
    [ValidateSet('load', 'clean', 'build', 'buildConda', 'debug', 'make', 'wfTest', 'wfTestConda', 'package', 'run')]
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
$RepositoryBuildRoot = Join-Path $RepoDir 'build'
$BuildRoot = Join-Path $RepositoryBuildRoot 'windows'
$ReleaseBuildDir = Join-Path $BuildRoot 'system'
$CondaBuildDir = Join-Path $BuildRoot 'conda'
$DebugBuildDir = Join-Path $BuildRoot 'conda-debug'
$PackageVersion = if ($env:SIPROS_VERSION) { $env:SIPROS_VERSION } else { '6.0.0' }
$PackageBaseName = "sipros_windows_$PackageVersion"
$InstallerProductName = 'Sipros'
$InstallerManufacturerName = 'Sipros'
$InstallerProductInfoUrl = 'https://github.com/xyz1396/sipros5'
$WixVersion = '5.0.2'

function Invoke-CheckedNative(
    [string]$Executable,
    [string[]]$ArgumentList,
    [string]$Description
) {
    & $Executable @ArgumentList
    $exitCode = $LASTEXITCODE
    if ($exitCode -ne 0) {
        throw "$Description failed with exit code $exitCode."
    }
}

function Wait-NinjaBuildDirectory([string]$Directory) {
    $ninjaLock = Join-Path $Directory '.ninja_lock'
    $deadline = [DateTime]::UtcNow.AddSeconds(60)
    $announced = $false
    while (Test-Path -LiteralPath $ninjaLock) {
        if (-not $announced) {
            Write-Warning "Another Ninja process is using $Directory; waiting for it to finish."
            $announced = $true
        }
        if ([DateTime]::UtcNow -ge $deadline) {
            throw "Another Ninja process is still using $Directory after 60 seconds. Stop the other build and retry."
        }
        Start-Sleep -Milliseconds 500
    }
}

function Invoke-CheckedCMakeBuild(
    [string]$Directory,
    [string[]]$ArgumentList,
    [string]$Description
) {
    for ($attempt = 1; $attempt -le 2; $attempt++) {
        Wait-NinjaBuildDirectory $Directory
        & cmake @ArgumentList
        $exitCode = $LASTEXITCODE
        if ($exitCode -eq 0) { return }

        $ninjaLock = Join-Path $Directory '.ninja_lock'
        if ($attempt -eq 1 -and (Test-Path -LiteralPath $ninjaLock)) {
            Write-Warning 'Ninja started concurrently after the lock check; waiting and retrying once.'
            continue
        }
        throw "$Description failed with exit code $exitCode."
    }
}

function Enter-WindowsBuildLock {
    $lockKey = $RepoDir.ToLowerInvariant() -replace '[^a-z0-9.-]', '_'
    $lockPath = Join-Path $env:TEMP ("sipros5-build-$lockKey.lock")
    try {
        return [System.IO.FileStream]::new(
            $lockPath,
            [System.IO.FileMode]::OpenOrCreate,
            [System.IO.FileAccess]::ReadWrite,
            [System.IO.FileShare]::None)
    } catch [System.IO.IOException] {
        throw "Another Sipros build, test, or package command is already running for $RepoDir. Wait for it to finish and retry."
    }
}

function Clear-WindowsOutputs {
    $binDirectory = Join-Path $RepoDir 'bin'
    if (Test-Path -LiteralPath $RepositoryBuildRoot) {
        $resolvedBuild = [System.IO.Path]::GetFullPath($RepositoryBuildRoot)
        $expectedBuild = [System.IO.Path]::GetFullPath((Join-Path $RepoDir 'build'))
        if ($resolvedBuild -ne $expectedBuild) {
            throw "Refusing to remove unexpected build directory: $resolvedBuild"
        }
        Remove-Item -LiteralPath $resolvedBuild -Recurse -Force
    }
    New-Item -ItemType Directory -Path $binDirectory -Force | Out-Null
    # Keep checked-in/downloaded runtime inputs (Raxport and DIA-NN models)
    # available for the next configure. Remove generated binaries from both
    # platform builds so either clean command leaves the same publish state.
    foreach ($path in @(
        (Join-Path $binDirectory 'aerith'),
        (Join-Path $binDirectory 'aerith.exe'),
        (Join-Path $binDirectory 'sipros'),
        (Join-Path $binDirectory 'sipros.exe'),
        (Join-Path $binDirectory 'siproswf'),
        (Join-Path $binDirectory 'siproswf.exe'),
        (Join-Path $binDirectory 'siprosMPI'),
        (Join-Path $binDirectory 'siprosMPI.exe')
    )) {
        if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path -Force }
    }
    $runtimeDirectory = [System.IO.Path]::GetFullPath((Join-Path $binDirectory 'lib'))
    $expectedRuntimeDirectory = [System.IO.Path]::GetFullPath(
        (Join-Path (Join-Path $RepoDir 'bin') 'lib'))
    if ($runtimeDirectory -ne $expectedRuntimeDirectory) {
        throw "Refusing to remove unexpected runtime directory: $runtimeDirectory"
    }
    if (Test-Path -LiteralPath $runtimeDirectory) {
        Remove-Item -LiteralPath $runtimeDirectory -Recurse -Force
    }
    foreach ($pattern in @(
        'siprosRelease.zip',
        'siprosRelease.msi',
        'sipros_windows_*.zip',
        'sipros_windows_*.msi',
        'sipros_linux_*.zip',
        'sipros_linux_*.AppImage'
    )) {
        foreach ($package in (Get-ChildItem -LiteralPath $RepoDir -File -Filter $pattern)) {
            if ([System.IO.Path]::GetFullPath($package.DirectoryName) -ne
                [System.IO.Path]::GetFullPath($RepoDir)) {
                throw "Refusing to remove release package outside the repository root: $($package.FullName)"
            }
            Remove-Item -LiteralPath $package.FullName -Force
        }
    }
}

function Get-ReleaseDotnetExecutable {
    $dotnetCandidates = @(
        (Join-Path $env:CONDA_PREFIX 'dotnet.exe'),
        (Join-Path $env:CONDA_PREFIX 'dotnet\dotnet.exe'),
        (Join-Path $env:CONDA_PREFIX 'Library\bin\dotnet.exe'),
        (Join-Path $env:CONDA_PREFIX 'Scripts\dotnet.exe')
    )
    $dotnet = $dotnetCandidates | Where-Object { Test-Path -LiteralPath $_ } | Select-Object -First 1
    if (-not $dotnet) {
        throw "The '$ReleaseEnvironmentName' environment must provide dotnet-sdk. Install it with: micromamba install -n $ReleaseEnvironmentName -c conda-forge 'dotnet-sdk=8.*'"
    }
    return $dotnet
}

function Get-WixExecutable {
    $dotnet = Get-ReleaseDotnetExecutable
    if (-not $env:CONDA_PREFIX) {
        throw "The '$ReleaseEnvironmentName' environment must be active before installing WiX."
    }
    $toolDirectory = Join-Path $env:CONDA_PREFIX "tools\wix-$WixVersion"
    $wix = Join-Path $toolDirectory 'wix.exe'
    if (-not (Test-Path -LiteralPath $wix)) {
        New-Item -ItemType Directory -Path $toolDirectory -Force | Out-Null
        Write-Host "Installing WiX $WixVersion in $toolDirectory with $dotnet"
        Invoke-CheckedNative -Executable $dotnet `
            -ArgumentList @(
                'tool', 'install', 'wix', '--version', $WixVersion,
                '--tool-path', $toolDirectory
            ) `
            -Description 'WiX tool installation' | Out-Host
        if (-not (Test-Path -LiteralPath $wix)) {
            throw "WiX installation did not create $wix"
        }
    }

    $uiExtension = Join-Path $toolDirectory ".wix\extensions\WixToolset.UI.wixext\$WixVersion\wixext5\WixToolset.UI.wixext.dll"
    if (-not (Test-Path -LiteralPath $uiExtension)) {
        Write-Host "Installing WiX UI extension $WixVersion in $toolDirectory"
        Push-Location -LiteralPath $toolDirectory
        try {
            Invoke-CheckedNative -Executable $wix `
                -ArgumentList @('extension', 'add', "WixToolset.UI.wixext/$WixVersion") `
                -Description 'WiX UI extension installation' | Out-Host
        } finally {
            Pop-Location
        }
        if (-not (Test-Path -LiteralPath $uiExtension)) {
            throw "WiX UI extension installation did not create $uiExtension"
        }
    }
    return $wix
}

function Assert-PackageVersion {
    $match = [regex]::Match($PackageVersion, '^(\d+)\.(\d+)\.(\d+)$')
    if (-not $match.Success) {
        throw "SIPROS_VERSION must be an MSI-compatible three-part numeric version, for example 6.0.0; got '$PackageVersion'."
    }
    if ([int]$match.Groups[1].Value -gt 255 -or
        [int]$match.Groups[2].Value -gt 255 -or
        [int]$match.Groups[3].Value -gt 65535) {
        throw "SIPROS_VERSION components exceed MSI limits (255.255.65535): '$PackageVersion'."
    }
}

function Convert-TextLicenseToRtf(
    [string]$Source,
    [string]$Destination
) {
    $text = [System.IO.File]::ReadAllText($Source, [System.Text.Encoding]::UTF8)
    $rtf = [System.Text.StringBuilder]::new()
    [void]$rtf.Append('{\rtf1\ansi\deff0{\fonttbl{\f0\fnil Segoe UI;}}\uc1\pard\f0\fs18 ')
    foreach ($character in $text.ToCharArray()) {
        switch ($character) {
            '\' { [void]$rtf.Append('\\') }
            '{' { [void]$rtf.Append('\{') }
            '}' { [void]$rtf.Append('\}') }
            "`r" { }
            "`n" { [void]$rtf.Append("\par`r`n") }
            "`t" { [void]$rtf.Append('\tab ') }
            default {
                $codePoint = [int]$character
                if ($codePoint -le 127) {
                    [void]$rtf.Append($character)
                } else {
                    if ($codePoint -gt 32767) { $codePoint -= 65536 }
                    [void]$rtf.Append("\u$codePoint?")
                }
            }
        }
    }
    [void]$rtf.Append('}')
    [System.IO.File]::WriteAllText(
        $Destination, $rtf.ToString(), [System.Text.Encoding]::ASCII)
}

function Import-VisualStudioEnvironment {
    if (Get-Command cl.exe -ErrorAction SilentlyContinue) { return }
    $vswhere = Join-Path ${env:ProgramFiles(x86)} 'Microsoft Visual Studio\Installer\vswhere.exe'
    if (-not (Test-Path -LiteralPath $vswhere)) {
        throw 'Visual Studio 2022 with the Desktop development with C++ workload is required.'
    }
    $installation = & $vswhere -latest -products '*' -requires Microsoft.VisualStudio.Component.VC.Tools.x86.x64 -property installationPath
    if (-not $installation) { throw 'Install the Visual Studio 2022 Desktop development with C++ workload.' }
    $devShellModule = Join-Path $installation 'Common7\Tools\Microsoft.VisualStudio.DevShell.dll'
    if (-not (Test-Path -LiteralPath $devShellModule)) {
        throw "Visual Studio's PowerShell developer-shell module was not found: $devShellModule"
    }
    Import-Module -Name $devShellModule
    Enter-VsDevShell -VsInstallPath $installation -SkipAutomaticLocation `
        -DevCmdArguments '-no_logo -arch=x64 -host_arch=x64'
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
    $check = "import glob,os,re,torch;p=os.environ['CONDA_PREFIX'];v=tuple(map(int,re.match(r'\d+(?:\.\d+)+',torch.__version__).group().split('.')));assert v$operator(2,12,1);assert glob.glob(os.path.join(p,'conda-meta','hdf5-2.*-nompi*.json'));assert glob.glob(os.path.join(p,'conda-meta','imgui-1.92.9-*.json'));assert glob.glob(os.path.join(p,'conda-meta','libvulkan-headers-*.json'))$cpuCheck$gpuCheck"
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
    $activation = & $MambaExe shell activate -n $Name -s powershell
    $exitCode = $LASTEXITCODE
    if ($exitCode -ne 0) {
        throw "Unable to activate micromamba environment '$Name' (exit code $exitCode)."
    }
    Invoke-Expression ($activation -join [Environment]::NewLine)
    if (-not (Test-ActiveMicromambaEnvironment $Name)) {
        throw "Micromamba activation did not select the expected '$Name' environment."
    }
    Write-Host "Activated micromamba environment '$Name'."
}

function Test-ActiveMicromambaEnvironment([string]$Name) {
    if (-not $env:CONDA_PREFIX) { return $false }
    $expectedPrefix = Join-Path $MambaRootPrefix "envs\$Name"
    try {
        return [System.IO.Path]::GetFullPath($env:CONDA_PREFIX).TrimEnd('\') -eq
            [System.IO.Path]::GetFullPath($expectedPrefix).TrimEnd('\')
    } catch {
        return $false
    }
}

function Invoke-WindowsBuild([string]$Directory, [string]$Configuration) {
    $torchPrefix = python -c 'import torch; print(torch.utils.cmake_prefix_path)'
    if ($LASTEXITCODE -ne 0) { throw 'Unable to locate the PyTorch CMake package.' }
    $torchDir = Join-Path $torchPrefix.Trim() 'Torch'
    $hdf5Dir = Join-Path $env:CONDA_PREFIX 'Library\cmake'
    $arguments = @(
        '-S', $RepoDir, '-B', $Directory, '-G', 'Ninja',
        "-DCMAKE_BUILD_TYPE=$Configuration", '-DBUILD_CONDA=ON',
        "-DSIPROS_VERSION=$PackageVersion",
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
        if (-not $env:CUDAARCHS) {
            $env:CUDAARCHS = '50-real;52-real;60-real;61-real;70-real;75-real;80-real;86-real;89-real;90-real;100-real;101-real;120'
        }
        if (-not $env:TORCH_CUDA_ARCH_LIST) {
            $env:TORCH_CUDA_ARCH_LIST = '5.0;5.2;6.0;6.1;7.0;7.5;8.0;8.6;8.9;9.0;10.0;10.1;12.0+PTX'
        }
        $cudaTargetInclude = Join-Path $env:CONDA_PREFIX 'Library\include\targets\x64'
        if ((Test-Path -LiteralPath $cudaTargetInclude) -and
            ($env:INCLUDE -split ';' -notcontains $cudaTargetInclude)) {
            $env:INCLUDE = "$env:INCLUDE;$cudaTargetInclude"
        }
        $arguments += '-DAERITH_ENABLE_CUDA=ON'
        $arguments += "-DCMAKE_CUDA_COMPILER=$($nvcc.Source)"
        $arguments += "-DCUDAToolkit_ROOT=$env:CONDA_PREFIX"
        Write-Host "Configuring CUDA with $($nvcc.Source)"
    } else {
        $arguments += '-DAERITH_ENABLE_CUDA=OFF'
    }
    Invoke-CheckedNative -Executable 'cmake' -ArgumentList $arguments `
        -Description 'CMake configuration'
    Invoke-CheckedCMakeBuild -Directory $Directory `
        -ArgumentList @(
            '--build', $Directory, '--target',
            'sipros', 'aerith', 'siproswf'
        ) `
        -Description 'Native build'
}

function Publish-WindowsBinaries([string]$SourceBinDirectory) {
    $executables = @('sipros.exe', 'aerith.exe', 'siproswf.exe')
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
    New-Item -ItemType Directory -Path $binDestination -Force | Out-Null
    foreach ($name in @($executables + $models)) {
        Copy-Item -LiteralPath (Join-Path $SourceBinDirectory $name) `
            -Destination (Join-Path $binDestination $name) -Force
    }
}

function Test-PackagedExecutable(
    [string]$Executable,
    [string[]]$ArgumentList
) {
    $startInfo = [System.Diagnostics.ProcessStartInfo]::new()
    $startInfo.FileName = $Executable
    $startInfo.WorkingDirectory = Split-Path -Parent $Executable
    $startInfo.UseShellExecute = $false
    $startInfo.CreateNoWindow = $true
    $startInfo.RedirectStandardOutput = $true
    $startInfo.RedirectStandardError = $true
    $startInfo.Environment['PATH'] = "$env:SystemRoot\System32;$env:SystemRoot"
    foreach ($argument in $ArgumentList) { $startInfo.ArgumentList.Add($argument) }

    $process = [System.Diagnostics.Process]::new()
    $process.StartInfo = $startInfo
    try {
        if (-not $process.Start()) { throw "Unable to start packaged executable: $Executable" }
        $stdoutTask = $process.StandardOutput.ReadToEndAsync()
        $stderrTask = $process.StandardError.ReadToEndAsync()
        if (-not $process.WaitForExit(30000)) {
            $process.Kill($true)
            throw "Packaged executable did not finish its startup check: $Executable"
        }
        $stdout = $stdoutTask.GetAwaiter().GetResult()
        $stderr = $stderrTask.GetAwaiter().GetResult()
        if ($process.ExitCode -ne 0) {
            $details = ($stdout + "`n" + $stderr).Trim()
            if ($details.Length -gt 2000) { $details = $details.Substring(0, 2000) }
            throw "Packaged executable failed its startup check with exit code $($process.ExitCode): $Executable`n$details"
        }
    } finally {
        $process.Dispose()
    }
}

function Copy-WindowsRuntimeDependencies(
    [string[]]$RootBinaries,
    [string]$Destination
) {
    $dumpbin = Get-Command dumpbin.exe -ErrorAction SilentlyContinue
    if (-not $dumpbin) {
        throw 'Visual Studio dumpbin.exe is required to resolve package DLL dependencies.'
    }

    $condaBin = Join-Path $env:CONDA_PREFIX 'Library\bin'
    if (-not (Test-Path -LiteralPath $condaBin)) {
        throw "Conda runtime directory was not found: $condaBin"
    }
    $condaDlls = @{}
    Get-ChildItem -LiteralPath $condaBin -Filter '*.dll' -File | ForEach-Object {
        $condaDlls[$_.Name] = $_.FullName
    }

    $queue = [System.Collections.Generic.Queue[string]]::new()
    $visited = [System.Collections.Generic.HashSet[string]]::new(
        [System.StringComparer]::OrdinalIgnoreCase)
    $runtimeDlls = [System.Collections.Generic.SortedDictionary[string,string]]::new(
        [System.StringComparer]::OrdinalIgnoreCase)
    foreach ($binary in $RootBinaries) { $queue.Enqueue($binary) }

    while ($queue.Count -gt 0) {
        $binary = $queue.Dequeue()
        if (-not $visited.Add($binary)) { continue }
        $imports = & $dumpbin.Source /nologo /dependents $binary 2>&1
        if ($LASTEXITCODE -ne 0) {
            throw "Unable to inspect DLL dependencies for $binary"
        }
        $dependencyNames = [System.Collections.Generic.HashSet[string]]::new(
            [System.StringComparer]::OrdinalIgnoreCase)
        foreach ($line in $imports) {
            if ($line -match '^\s+([^\s]+\.dll)\s*$') {
                [void]$dependencyNames.Add($Matches[1])
            }
        }
        # Conda's SDL2 compatibility library loads SDL3 at runtime instead of
        # listing it in the PE import table.
        if ([System.IO.Path]::GetFileName($binary) -ieq 'SDL2.dll') {
            [void]$dependencyNames.Add('SDL3.dll')
        }
        # oneMKL selects both a general CPU kernel and a VML kernel at runtime.
        # Include every shipped ISA alternative so the MSI works on supported
        # CPUs without relying on DLLs from an existing Conda installation.
        if ([System.IO.Path]::GetFileName($binary) -ieq 'mkl_core.3.dll') {
            foreach ($name in @(
                'mkl_avx10.3.dll',
                'mkl_avx2.3.dll',
                'mkl_avx512.3.dll',
                'mkl_def.3.dll',
                'mkl_mc3.3.dll',
                'mkl_vml_avx10.3.dll',
                'mkl_vml_avx2.3.dll',
                'mkl_vml_avx512.3.dll',
                'mkl_vml_cmpt.3.dll',
                'mkl_vml_def.3.dll',
                'mkl_vml_mc3.3.dll'
            )) {
                [void]$dependencyNames.Add($name)
            }
        }
        foreach ($name in $dependencyNames) {
            $adjacentDll = Join-Path (Split-Path -Parent $binary) $name
            $source = if (Test-Path -LiteralPath $adjacentDll) {
                $adjacentDll
            } elseif ($condaDlls.ContainsKey($name)) {
                $condaDlls[$name]
            } else {
                $null
            }
            if ($source) {
                if (-not $runtimeDlls.ContainsKey($name)) {
                    $runtimeDlls.Add($name, $source)
                    $queue.Enqueue($source)
                }
                continue
            }
            $systemDll = Join-Path $env:SystemRoot "System32\$name"
            if ((Test-Path -LiteralPath $systemDll) -or $name -match '^(api|ext)-ms-win-') {
                continue
            }
            throw "Unresolved non-system DLL '$name' imported by $binary"
        }
    }

    foreach ($entry in $runtimeDlls.GetEnumerator()) {
        Copy-Item -LiteralPath $entry.Value -Destination (Join-Path $Destination $entry.Key)
    }
    Write-Host "Packaged $($runtimeDlls.Count) linked and runtime-selected DLLs."
}

function New-WindowsPackage([string]$SourceBinDirectory) {
    Assert-PackageVersion
    $stageRoot = Join-Path $env:TEMP ("sipros5-package-" + [guid]::NewGuid().ToString('N'))
    $stage = Join-Path $stageRoot 'sipros'
    $archive = Join-Path $RepoDir "$PackageBaseName.zip"
    $installer = Join-Path $RepoDir "$PackageBaseName.msi"
    try {
        $lib = Join-Path $stage 'lib'
        New-Item -ItemType Directory -Path $lib -Force | Out-Null
        Copy-Item -LiteralPath (Join-Path $RepoDir 'LICENSE') -Destination $lib
        $workflow = Join-Path $SourceBinDirectory 'siproswf.exe'
        if (-not (Test-Path -LiteralPath $workflow)) { throw "Missing built runtime file: $workflow" }
        Copy-Item -LiteralPath $workflow -Destination $stage
        foreach ($name in 'sipros.exe', 'aerith.exe') {
            $source = Join-Path $SourceBinDirectory $name
            if (-not (Test-Path -LiteralPath $source)) { throw "Missing built runtime file: $source" }
            Copy-Item -LiteralPath $source -Destination $lib
        }
        $raxport = Join-Path $RepoDir 'bin\Raxport-win-x64.exe'
        if (-not (Test-Path -LiteralPath $raxport)) { throw "Missing runtime binary: $raxport" }
        Copy-Item -LiteralPath $raxport -Destination $lib
        foreach ($name in 'diann-2.6.1-fragmentation.pt', 'diann-2.6.1-retention-time.pt') {
            $source = Join-Path $SourceBinDirectory $name
            if (-not (Test-Path -LiteralPath $source)) { throw "Missing model file: $source" }
            Copy-Item -LiteralPath $source -Destination $lib
        }
        Copy-WindowsRuntimeDependencies `
            -RootBinaries @(
                $workflow,
                (Join-Path $SourceBinDirectory 'sipros.exe'),
                (Join-Path $SourceBinDirectory 'aerith.exe'),
                $raxport
            ) `
            -Destination $lib
        foreach ($name in @(
            'hdf5.dll', 'hdf5_cpp.dll', 'torch_cpu.dll', 'c10.dll',
            'libomp.dll', 'imgui.dll', 'glfw3.dll'
        )) {
            $runtime = Join-Path $lib $name
            if (-not (Test-Path -LiteralPath $runtime)) {
                throw "Missing required packaged runtime DLL: $runtime"
            }
        }
        Test-PackagedExecutable (Join-Path $lib 'sipros.exe') @('--help')
        Test-PackagedExecutable (Join-Path $lib 'aerith.exe') @('--help')
        Test-PackagedExecutable (Join-Path $stage 'siproswf.exe') @('--help')
        Write-Host 'Packaged Sipros, Aerith, and Sipros workflow runtime checks passed.'
        if (Test-Path -LiteralPath $archive) { Remove-Item -LiteralPath $archive -Force }
        Compress-Archive -LiteralPath $stage -DestinationPath $archive -CompressionLevel Optimal
        Write-Host "Package created: $archive"

        $wix = Get-WixExecutable
        $wixUiExtension = Join-Path (Split-Path -Parent $wix) ".wix\extensions\WixToolset.UI.wixext\$WixVersion\wixext5\WixToolset.UI.wixext.dll"
        $licenseDirectory = Join-Path $stageRoot 'license-ui'
        New-Item -ItemType Directory -Path $licenseDirectory -Force | Out-Null
        Convert-TextLicenseToRtf `
            -Source (Join-Path $RepoDir 'LICENSE') `
            -Destination (Join-Path $licenseDirectory 'LICENSE.rtf')
        $wixIntermediate = Join-Path $stageRoot 'wix'
        New-Item -ItemType Directory -Path $wixIntermediate -Force | Out-Null
        if (Test-Path -LiteralPath $installer) { Remove-Item -LiteralPath $installer -Force }
        Invoke-CheckedNative -Executable $wix `
            -ArgumentList @(
                'build', (Join-Path $RepoDir 'wf33\sipros.wxs'),
                '-arch', 'x64',
                '-bindpath', "Stage=$stage",
                '-bindpath', "Source=$RepoDir",
                '-bindpath', "License=$licenseDirectory",
                '-ext', $wixUiExtension,
                '-define', "PackageVersion=$PackageVersion",
                '-define', "ProductName=$InstallerProductName",
                '-define', "ManufacturerName=$InstallerManufacturerName",
                '-define', "ProductInfoUrl=$InstallerProductInfoUrl",
                '-intermediateFolder', $wixIntermediate,
                '-pdbtype', 'none',
                '-out', $installer
            ) `
            -Description 'MSI package build'
        Write-Host "Installer created: $installer"
    } finally {
        if (Test-Path -LiteralPath $stageRoot) { Remove-Item -LiteralPath $stageRoot -Recurse -Force }
    }
}

$env:MAMBA_ROOT_PREFIX = $MambaRootPrefix
if ($Command -eq 'clean') {
    $buildLock = Enter-WindowsBuildLock
    try { Clear-WindowsOutputs } finally { $buildLock.Dispose() }
    return
}
if ($Command -eq 'load') { return }

Import-VisualStudioEnvironment
if (-not $EnvironmentActive) {
    if ($Command -in 'build', 'package', 'make', 'wfTest', 'run') {
        Require-ReleaseEnvironment
        if (Test-ActiveMicromambaEnvironment $ReleaseEnvironmentName) {
            Write-Host "Using active micromamba environment '$ReleaseEnvironmentName'."
        } else {
            Enter-MicromambaEnvironment $ReleaseEnvironmentName
        }
    } elseif ($Command -in 'buildConda', 'debug', 'wfTestConda') {
        Require-CondaEnvironment
        if (Test-ActiveMicromambaEnvironment $CondaEnvironmentName) {
            Write-Host "Using active micromamba environment '$CondaEnvironmentName'."
        } else {
            Enter-MicromambaEnvironment $CondaEnvironmentName
        }
    }
}

$buildLock = Enter-WindowsBuildLock
try {
switch ($Command) {
    'build' {
        Invoke-WindowsBuild $ReleaseBuildDir 'Release'
        Publish-WindowsBinaries (Join-Path $ReleaseBuildDir 'bin')
        Write-Host "Windows release binaries published in $RepoDir\bin"
    }
    'buildConda' {
        Invoke-WindowsBuild $CondaBuildDir 'Release'
        Publish-WindowsBinaries (Join-Path $CondaBuildDir 'bin')
        Write-Host "Windows Conda binaries published in $RepoDir\bin"
    }
    'debug' { Invoke-WindowsBuild $DebugBuildDir 'Debug' }
    'make' {
        Invoke-CheckedCMakeBuild -Directory $ReleaseBuildDir `
            -ArgumentList @('--build', $ReleaseBuildDir) `
            -Description 'Native build'
    }
    'wfTest' {
        Invoke-CheckedCMakeBuild -Directory $ReleaseBuildDir `
            -ArgumentList @('--build', $ReleaseBuildDir, '--target', 'siproswf_test') `
            -Description 'siproswf test build'
        Invoke-CheckedNative -Executable 'ctest' `
            -ArgumentList @('--test-dir', $ReleaseBuildDir, '-R', '^siproswf_core$', '--output-on-failure') `
            -Description 'siproswf tests'
    }
    'wfTestConda' {
        Invoke-CheckedCMakeBuild -Directory $CondaBuildDir `
            -ArgumentList @('--build', $CondaBuildDir, '--target', 'siproswf_test') `
            -Description 'siproswf Conda test build'
        Invoke-CheckedNative -Executable 'ctest' `
            -ArgumentList @('--test-dir', $CondaBuildDir, '-R', '^siproswf_core$', '--output-on-failure') `
            -Description 'siproswf Conda tests'
    }
    'package' {
        [void](Get-ReleaseDotnetExecutable)
        Invoke-WindowsBuild $ReleaseBuildDir 'Release'
        Publish-WindowsBinaries (Join-Path $ReleaseBuildDir 'bin')
        New-WindowsPackage (Join-Path $ReleaseBuildDir 'bin')
    }
    'run' {
        Write-Host 'Open the ImGui workflow or run a headless search, for example:'
        Write-Host '  micromamba run -n sipros5 bin\siproswf.exe'
        Write-Host '  micromamba run -n sipros5 bin\siproswf.exe --regular-fasta-search -i input.h5 -f db.faa -o output'
    }
}
} finally {
    $buildLock.Dispose()
}
