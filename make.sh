#!/bin/bash
set -e

# Build dependency setup
#
# `build` and `package` run in a dedicated CPU release environment. The
# environment is created automatically when it is missing or has stale pins:
#   micromamba create -n sipros5-release -c conda-forge \
#     sysroot_linux-64=2.17 gcc_linux-64 gxx_linux-64 cmake ninja patchelf \
#     "hdf5=2.*=nompi*" \
#     "pytorch-cpu=2.12.1=cpu_mkl*"
# 
# Micromamba/Conda build (./make.sh buildConda):
#   micromamba create -n sipros5 -c conda-forge \
#     hdf5 h5py openmpi cmake ninja gcc_linux-64 gxx_linux-64 gdb gperftools \
#     python=3.12 lxml pandas sysroot_linux-64=2.34 scikit-learn matplotlib \
#     pytorch-gpu=2.12.1 cuda-cudart-dev=12.9 libcublas-dev=12.9 \
#     cuda-nvrtc-dev=12.9 cuda-nvcc=12.9 cuda-nvtx-dev=12.9
#
# Activate the environment before running Conda-built binaries so their
# OpenMP, compiler-runtime, and MPI shared libraries are available:
#   micromamba activate sipros5
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAMBA_EXE="${MAMBA_EXE:-micromamba}"
if ! command -v "$MAMBA_EXE" >/dev/null 2>&1 && [ -x "$HOME/.local/bin/micromamba" ]; then
    MAMBA_EXE="$HOME/.local/bin/micromamba"
fi
RELEASE_ENV_NAME="${RELEASE_ENV_NAME:-sipros5-release}"
RELEASE_ENV_ACTIVE="${SIPROS_RELEASE_ENV_ACTIVE:-0}"

ensure_release_environment() {
    if "$MAMBA_EXE" run -n "$RELEASE_ENV_NAME" python -c \
        'import glob, os, torch
prefix = os.environ["CONDA_PREFIX"]
assert torch.__version__.split("+")[0] == "2.12.1"
assert glob.glob(prefix + "/conda-meta/pytorch-2.12.1-cpu_mkl*.json")
assert glob.glob(prefix + "/conda-meta/sysroot_linux-64-2.17-*.json")
assert glob.glob(prefix + "/conda-meta/gcc_linux-64-*.json")
assert glob.glob(prefix + "/conda-meta/gxx_linux-64-*.json")
assert glob.glob(prefix + "/conda-meta/hdf5-2.*-nompi*.json")
assert os.path.isfile(prefix + "/lib/libhdf5.so")
assert os.path.isfile(prefix + "/lib/libhdf5_cpp.so")
assert os.access(prefix + "/bin/patchelf", os.X_OK)' \
        >/dev/null 2>&1; then
        return
    fi

    local action=create
    if "$MAMBA_EXE" run -n "$RELEASE_ENV_NAME" true >/dev/null 2>&1; then
        action=install
    fi
    "$MAMBA_EXE" "$action" -y -n "$RELEASE_ENV_NAME" \
        -c conda-forge --strict-channel-priority \
        "sysroot_linux-64=2.17" \
        gcc_linux-64 gxx_linux-64 cmake ninja patchelf \
        "hdf5=2.*=nompi*" \
        "pytorch-cpu=2.12.1=cpu_mkl*"
}

# Only the release build and package commands are isolated in the pinned
# environment. All other commands retain their existing execution paths.
if [[ "${1:-}" = "build" || "${1:-}" = "package" ]] &&
   [ "$RELEASE_ENV_ACTIVE" != "1" ]; then
    ensure_release_environment
    exec "$MAMBA_EXE" run -n "$RELEASE_ENV_NAME" \
        env -u LD_LIBRARY_PATH -u LDFLAGS SIPROS_RELEASE_ENV_ACTIVE=1 \
        "$REPO_DIR/make.sh" "$@"
fi

if [ "$RELEASE_ENV_ACTIVE" = "1" ]; then
    RELEASE_PREFIX="${CONDA_PREFIX:?Missing release environment prefix}"
    SYSTEM_CMAKE="$RELEASE_PREFIX/bin/cmake"
    SYSTEM_NINJA="$RELEASE_PREFIX/bin/ninja"
    SYSTEM_CC="$RELEASE_PREFIX/bin/x86_64-conda-linux-gnu-gcc"
    SYSTEM_CXX="$RELEASE_PREFIX/bin/x86_64-conda-linux-gnu-g++"
    SYSTEM_TORCH_ROOT="$RELEASE_PREFIX"
    if [ -f "$RELEASE_PREFIX/lib/cmake/hdf5/hdf5-config.cmake" ]; then
        RELEASE_HDF5_DIR="$RELEASE_PREFIX/lib/cmake/hdf5"
    else
        # Older Conda HDF5 packages install their CMake package here.
        RELEASE_HDF5_DIR="$RELEASE_PREFIX/cmake"
    fi
else
    SYSTEM_CMAKE="${SYSTEM_CMAKE:-/usr/bin/cmake}"
    SYSTEM_NINJA="${SYSTEM_NINJA:-/usr/bin/ninja}"
    SYSTEM_CC="${SYSTEM_CC:-/usr/bin/gcc}"
    SYSTEM_CXX="${SYSTEM_CXX:-/usr/bin/g++}"
    SYSTEM_TORCH_ROOT="${SYSTEM_TORCH_ROOT:-/opt/libtorch}"
fi
SYSTEM_TORCH_LIB_DIR="$SYSTEM_TORCH_ROOT/lib"
BUILD_ROOT="$REPO_DIR/build"
SYSTEM_BUILD_DIR="$BUILD_ROOT/system"
CONDA_BUILD_DIR="$BUILD_ROOT/conda"
CONDA_DEBUG_BUILD_DIR="$BUILD_ROOT/conda-debug"
DIANN_MODEL_NAME="diann-2.6.1-fragmentation.pt"
DIANN_RT_MODEL_NAME="diann-2.6.1-retention-time.pt"
AERITH_MIN_TORCH_VERSION="2.12.1"
RUNTIME_TOOL_BINARIES=(
    Raxport-linux-x64
    sipros
    siprosMPI
    aerith
)
MKL_DISPATCH_LIBRARIES=(
    libmkl_avx2.so.3
    libmkl_avx512.so.3
    libmkl_def.so.3
    libmkl_vml_avx2.so.3
    libmkl_vml_avx512.so.3
    libmkl_vml_def.so.3
)

# Archive extraction and file copies can drop execute bits from the bundled
# tools.  Normalize every runtime binary before any build/run action so that
# `./make.sh load` can also repair an existing checkout without rebuilding it.
ensure_runtime_tool_permissions() {
    local binary path
    for binary in "${RUNTIME_TOOL_BINARIES[@]}"; do
        path="$REPO_DIR/tools/$binary"
        if [ -f "$path" ]; then
            chmod 0755 "$path"
        fi
    done
}

require_executable() {
    if [ ! -x "$1" ]; then
        echo "Missing required executable: $1" >&2
        exit 1
    fi
}

prepare_release_build_dir() {
    mkdir -p "$1"
    local cache="$1/CMakeCache.txt"
    local hdf5_dir="=$RELEASE_HDF5_DIR"
    local compiler="=$SYSTEM_CXX"
    local torch_root="=$SYSTEM_TORCH_ROOT"
    if [ -f "$cache" ] &&
       { ! grep -F 'HDF5_DIR:' "$cache" | grep -Fq "$hdf5_dir" ||
         ! grep -F 'CMAKE_CXX_COMPILER:' "$cache" | grep -Fq "$compiler" ||
         ! grep -F 'AERITH_TORCH_ROOT:' "$cache" | grep -Fq "$torch_root" ||
         ! grep -F 'AERITH_TORCH_CPU_ONLY:' "$cache" | grep -Fq '=ON'; }; then
        echo "Reconfiguring $1 with the dynamic release environment"
        rm -rf "$cache" "$1/CMakeFiles"
    fi
}

require_release_libraries() {
    local required
    if [ ! -f "$SYSTEM_TORCH_ROOT/include/torch/script.h" ]; then
        echo "Missing release LibTorch headers: $SYSTEM_TORCH_ROOT/include" >&2
        exit 1
    fi
    for required in libtorch_cpu.so libc10.so; do
        if [ ! -f "$SYSTEM_TORCH_LIB_DIR/$required" ]; then
            echo "Missing release CPU LibTorch library: $SYSTEM_TORCH_LIB_DIR/$required" >&2
            exit 1
        fi
    done
    for required in libhdf5.so libhdf5_cpp.so libgomp.so libstdc++.so.6; do
        if [ ! -f "$SYSTEM_TORCH_LIB_DIR/$required" ]; then
            echo "Missing dynamic release library: $SYSTEM_TORCH_LIB_DIR/$required" >&2
            exit 1
        fi
    done
    if [ ! -f "$RELEASE_HDF5_DIR/hdf5-config.cmake" ]; then
        echo "Missing release HDF5 CMake package: $RELEASE_HDF5_DIR" >&2
        exit 1
    fi
}

prepare_conda_build_dir() {
    mkdir -p "$1"
    local cache="$1/CMakeCache.txt"
    local hdf5_dir="=$CONDA_PREFIX/cmake"
    local torch_dir="=$TORCH_CMAKE_DIR"
    local torch_version="=$TORCH_VERSION"
    local cuda_compiler="=$CONDA_PREFIX/bin/nvcc"
    if [ -f "$cache" ] &&
       { ! grep -F 'HDF5_DIR:' "$cache" | grep -Fq "$hdf5_dir" ||
         ! grep -F 'Torch_DIR:' "$cache" | grep -Fq "$torch_dir" ||
         ! grep -F 'AERITH_TORCH_PACKAGE_VERSION:' "$cache" | grep -Fq "$torch_version" ||
         { [ -x "$CONDA_PREFIX/bin/nvcc" ] &&
           ! grep -F 'CMAKE_CUDA_COMPILER:' "$cache" | grep -Fq "$cuda_compiler"; }; }; then
        echo "Reconfiguring $1 with the current Conda HDF5, CUDA, and LibTorch"
        rm -rf "$cache" "$1/CMakeFiles"
    fi
}

run_sipros5() {
    "$MAMBA_EXE" run -n sipros5 "$@"
}

resolve_torch_cmake() {
    local torch_info torch_prefix
    torch_info=$(run_sipros5 python -c \
        'import re, sys, torch
required = tuple(map(int, sys.argv[1].split(".")))
installed = tuple(map(int, re.match(r"\d+(?:\.\d+)+", torch.__version__).group().split(".")))
installed += (0,) * (len(required) - len(installed))
if installed < required:
    raise SystemExit(f"PyTorch {sys.argv[1]} or newer is required; found {torch.__version__}")
print(torch.__version__.split("+")[0] + "\t" + torch.utils.cmake_prefix_path)' \
        "$AERITH_MIN_TORCH_VERSION") || {
        echo "PyTorch $AERITH_MIN_TORCH_VERSION or newer is required in the sipros5 Conda environment." >&2
        exit 1
    }
    TORCH_VERSION="${torch_info%%$'\t'*}"
    torch_prefix="${torch_info#*$'\t'}"
    TORCH_CMAKE_DIR="$torch_prefix/Torch"
    if [ ! -f "$TORCH_CMAKE_DIR/TorchConfig.cmake" ]; then
        echo "Missing dynamic LibTorch CMake package: $TORCH_CMAKE_DIR" >&2
        exit 1
    fi
}

# Ignore broken MPI wrappers left behind when a Conda environment is moved.
select_mpi_cxx() {
    local dir candidate
    local old_ifs="$IFS"
    IFS=:
    for dir in $PATH; do
        [ -n "$dir" ] || dir=.
        candidate="$dir/mpicxx"
        if [ -x "$candidate" ] && "$candidate" --showme:link >/dev/null 2>&1; then
            IFS="$old_ifs"
            printf '%s\n' "$candidate"
            return 0
        fi
    done
    IFS="$old_ifs"
    return 1
}

configure_mpi() {
    local mpi_cxx
    case " ${cmake_args[*]} " in
        *" -DMPI_CXX_COMPILER="*) return 0 ;;
    esac
    mpi_cxx=$(select_mpi_cxx) || {
        echo "No working MPI C++ wrapper (mpicxx) found on PATH." >&2
        exit 1
    }
    cmake_args+=("-DMPI_CXX_COMPILER=$mpi_cxx")
    # Avoid requiring legacy MPI-2 C++ bindings during CMake's probe.
    cmake_args+=("-DMPI_CXX_SKIP_MPICXX=ON")
}

verify_fully_dynamic_conda() {
    local binary="$1/sipros"
    local mpi_binary="$1/siprosMPI"
    local dependencies mpi_dependencies
    dependencies=$(ldd "$binary" 2>&1) || {
        echo "Expected a dynamically linked Conda OpenMP binary: $binary" >&2
        echo "$dependencies" >&2
        return 1
    }
    if ! grep -Eq 'libgomp\.so' <<<"$dependencies"; then
        echo "Conda OpenMP binary is not dynamically linked to libgomp: $binary" >&2
        echo "$dependencies" >&2
        return 1
    fi
    if ! grep -Eq 'libhdf5(_cpp)?\.so' <<<"$dependencies"; then
        echo "Conda OpenMP binary is not dynamically linked to HDF5: $binary" >&2
        echo "$dependencies" >&2
        return 1
    fi
    mpi_dependencies=$(ldd "$mpi_binary" 2>&1) || {
        echo "Expected a dynamically linked Conda MPI binary: $mpi_binary" >&2
        echo "$mpi_dependencies" >&2
        return 1
    }
    if ! grep -Eq 'libmpi\.so' <<<"$mpi_dependencies"; then
        echo "Conda MPI binary is not dynamically linked to MPI: $mpi_binary" >&2
        echo "$mpi_dependencies" >&2
        return 1
    fi
}

verify_dynamic_torch() {
    local binary="$1"
    local dependencies
    dependencies=$(ldd "$binary" 2>&1) || {
        echo "Expected Aerith to dynamically link LibTorch: $binary" >&2
        echo "$dependencies" >&2
        return 1
    }
    if grep -Fq 'not found' <<<"$dependencies" ||
       ! grep -Eq 'libtorch(_cpu)?\.so' <<<"$dependencies" ||
       ! grep -Eq 'libc10\.so' <<<"$dependencies"; then
        echo "Aerith does not have a complete dynamic LibTorch linkage: $binary" >&2
        echo "$dependencies" >&2
        return 1
    fi
}

verify_release_dynamic_binary() {
    local binary="$1"
    shift
    local dependencies required
    dependencies=$(LD_LIBRARY_PATH="$SYSTEM_TORCH_LIB_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
        ldd "$binary" 2>&1) || {
        echo "Expected a dynamically linked release binary: $binary" >&2
        echo "$dependencies" >&2
        return 1
    }
    if grep -Fq 'not found' <<<"$dependencies" ||
       grep -Eqi 'lib(torch_cuda|c10_cuda|cu(da|dnn|blas|pti|sparse|fft|rand|solver))' \
           <<<"$dependencies"; then
        echo "Release binary has an incomplete or non-CPU runtime: $binary" >&2
        echo "$dependencies" >&2
        return 1
    fi
    for required in "$@"; do
        if ! grep -Eq "$required" <<<"$dependencies"; then
            echo "Release binary is missing dynamic dependency $required: $binary" >&2
            echo "$dependencies" >&2
            return 1
        fi
    done
}

verify_glibc_217() {
    local binary="$1"
    local library_dir="${2:-}"
    local current highest major minor
    local -a files=("$binary")
    if [ -n "$library_dir" ]; then
        while IFS= read -r -d '' current; do
            files+=("$current")
        done < <(find "$library_dir" -maxdepth 1 -type f -print0)
    fi

    for current in "${files[@]}"; do
        highest=$(readelf --version-info "$current" 2>/dev/null |
            grep -oE 'GLIBC_[0-9]+\.[0-9]+' |
            sed 's/^GLIBC_//' | sort -Vu | tail -1)
        [ -n "$highest" ] || continue
        major="${highest%%.*}"
        minor="${highest#*.}"
        if (( major > 2 || (major == 2 && minor > 17) )); then
            echo "$current requires GLIBC_$highest; release limit is GLIBC_2.17" >&2
            return 1
        fi
    done
}

verify_packaged_release_binary() {
    local binary="$1"
    local library_dir="$2"
    local require_mkl="${3:-0}"
    local dependencies dynamic dependency arrow resolved remainder required
    dependencies=$(env -u LD_LIBRARY_PATH LD_LIBRARY_PATH="$library_dir" \
        ldd "$binary" 2>&1) || {
        echo "Unable to resolve packaged release runtime: $binary" >&2
        echo "$dependencies" >&2
        return 1
    }
    dynamic=$(readelf -d "$binary")
    if grep -Fq 'not found' <<<"$dependencies" ||
       grep -Eqi 'lib(torch_cuda|c10_cuda|cu(da|dnn|blas|pti|sparse|fft|rand|solver))' \
           <<<"$dependencies" ||
       grep -Fq "$SYSTEM_TORCH_LIB_DIR" <<<"$dynamic" ||
       ! grep -Fq '(RPATH)' <<<"$dynamic" ||
       ! grep -Fq '$ORIGIN/lib' <<<"$dynamic"; then
        echo "Packaged release runtime is incomplete or not relocatable: $binary" >&2
        echo "$dependencies" >&2
        return 1
    fi

    while read -r dependency arrow resolved remainder; do
        [ "$arrow" = "=>" ] || continue
        case "$resolved" in
            "$library_dir"/*) ;;
            /lib/*|/lib64/*|/usr/lib/*|/usr/lib64/*)
                case "$dependency" in
                    libc.so.*|libm.so.*|libpthread.so.*|libdl.so.*|librt.so.*|\
                    libresolv.so.*|libutil.so.*|libnsl.so.*|libanl.so.*) ;;
                    *)
                        echo "Packaged runtime uses host library $dependency: $resolved" >&2
                        return 1
                        ;;
                esac
                ;;
        esac
    done <<<"$dependencies"

    if [ "$require_mkl" = "1" ]; then
        for required in "${MKL_DISPATCH_LIBRARIES[@]}"; do
            if [ ! -f "$library_dir/$required" ]; then
                echo "Packaged Aerith is missing oneMKL dispatch library: $required" >&2
                return 1
            fi
        done
    fi
}

stage_release_runtime() {
    local binary="$1"
    local destination="$2"
    local current dependency resolved name needed source
    local -a queue=("$binary")
    local -A visited=()

    mkdir -p "$destination"
    while [ "${#queue[@]}" -gt 0 ]; do
        current="${queue[0]}"
        queue=("${queue[@]:1}")
        if [ -n "${visited[$current]+x}" ]; then
            continue
        fi
        visited["$current"]=1

        # ldd suppresses a second name when two SONAMEs resolve to the same
        # file (Conda's libgomp.so.1 -> libomp.so is one example). Read direct
        # NEEDED entries as well so every runtime name is present in the
        # relocatable package.
        while read -r needed; do
            source="$SYSTEM_TORCH_LIB_DIR/$needed"
            if [ -f "$source" ] && [ ! -f "$destination/$needed" ]; then
                if [[ "$needed" =~ (torch_cuda|c10_cuda|cuda|cudnn|cublas|cupti|cusparse|cufft|curand|cusolver) ]]; then
                    echo "CPU-only package unexpectedly requires $needed" >&2
                    return 1
                fi
                cp -L "$source" "$destination/$needed"
                chmod 0755 "$destination/$needed"
                queue+=("$source")
            fi
        done < <(readelf -d "$current" 2>/dev/null |
            sed -n 's/.*Shared library: \[\([^]]*\)\].*/\1/p')

        while read -r dependency resolved; do
            if [ "$resolved" = "not" ]; then
                echo "Unresolved runtime dependency while packaging $current: $dependency" >&2
                return 1
            fi
            case "$resolved" in
                "$SYSTEM_TORCH_LIB_DIR"/*)
                    name=$(basename "$resolved")
                    if [[ "$name" =~ (torch_cuda|c10_cuda|cuda|cudnn|cublas|cupti|cusparse|cufft|curand|cusolver) ]]; then
                        echo "CPU-only package unexpectedly requires $name" >&2
                        return 1
                    fi
                    if [ ! -f "$destination/$name" ]; then
                        cp -L "$resolved" "$destination/$name"
                        chmod 0755 "$destination/$name"
                    fi
                    queue+=("$resolved")
                    ;;
            esac
        done < <(LD_LIBRARY_PATH="$SYSTEM_TORCH_LIB_DIR${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
            ldd "$current" | awk '/=>/ { print $1, $3 }')
    done

    # oneMKL selects these CPU kernels with dlopen, so they do not appear in
    # ELF NEEDED records or ldd output. Bundle optimized kernels for the
    # prevalent AVX2/AVX-512 nodes plus the generic fallback.
    for name in "${MKL_DISPATCH_LIBRARIES[@]}"; do
        source="$SYSTEM_TORCH_LIB_DIR/$name"
        if [ ! -f "$source" ]; then
            echo "Missing required oneMKL dispatch library: $source" >&2
            return 1
        fi
        if [ ! -f "$destination/$name" ]; then
            cp -L "$source" "$destination/$name"
            chmod 0755 "$destination/$name"
        fi
    done
}

stage_publish_tools() {
    local source_bin_dir="$1"
    local profile="${2:-full}"
    local binaries
    if [ "$profile" = "cpu-release" ]; then
        binaries=(sipros aerith)
    else
        binaries=(sipros siprosMPI aerith)
    fi
    local assets=("$DIANN_MODEL_NAME" "$DIANN_RT_MODEL_NAME")
    local destinations=("$REPO_DIR/bin" "$REPO_DIR/tools")
    local asset binary destination tmpdir

    # Validate the complete build before replacing any published executable.
    for binary in "${binaries[@]}"; do
        if [ ! -x "$source_bin_dir/$binary" ]; then
            echo "Missing built executable: $source_bin_dir/$binary" >&2
            return 1
        fi
    done
    for asset in "${assets[@]}"; do
        if [ ! -f "$source_bin_dir/$asset" ]; then
            echo "Missing built runtime asset: $source_bin_dir/$asset" >&2
            return 1
        fi
    done

    for destination in "${destinations[@]}"; do
        mkdir -p "$destination"
        tmpdir=$(mktemp -d "$destination/.stage.XXXXXX")
        for binary in "${binaries[@]}"; do
            install -m 0755 "$source_bin_dir/$binary" "$tmpdir/$binary"
        done
        for asset in "${assets[@]}"; do
            install -m 0644 "$source_bin_dir/$asset" "$tmpdir/$asset"
        done
        for binary in "${binaries[@]}"; do
            mv -f "$tmpdir/$binary" "$destination/$binary"
        done
        if [ "$profile" = "cpu-release" ]; then
            rm -f "$destination/siprosMPI"
        fi
        for asset in "${assets[@]}"; do
            mv -f "$tmpdir/$asset" "$destination/$asset"
        done
        rmdir "$tmpdir"
        if [ "$profile" = "cpu-release" ]; then
            rm -rf "$destination/lib"
            stage_release_runtime "$destination/sipros" "$destination/lib"
            stage_release_runtime "$destination/aerith" "$destination/lib"
        fi
    done
    ensure_runtime_tool_permissions
}

cmake_args=()
if [ -n "${CMAKE_ARGS:-}" ]; then
    cmake_args=(${CMAKE_ARGS})
fi

# Repair bundled tool permissions even when no compilation is requested.
ensure_runtime_tool_permissions

case $1 in
"load") ;;
"clean")
    rm -rf "$BUILD_ROOT"
    mkdir "$BUILD_ROOT"
    rm -rf bin
    mkdir bin
    rm -f \
        "$REPO_DIR/tools/aerith" \
        "$REPO_DIR/tools/sipros" \
        "$REPO_DIR/tools/siprosMPI" \
        "$REPO_DIR/siprosRelease.zip"
    rm -rf "$REPO_DIR/tools/lib"
    ;;
"build")
    require_executable "$SYSTEM_CMAKE"
    require_executable "$SYSTEM_NINJA"
    require_executable "$SYSTEM_CC"
    require_executable "$SYSTEM_CXX"
    require_release_libraries
    prepare_release_build_dir "$SYSTEM_BUILD_DIR"
    mkdir -p tools
    cd "$SYSTEM_BUILD_DIR"
    CC="$SYSTEM_CC" CXX="$SYSTEM_CXX" "$SYSTEM_CMAKE" -G Ninja "${cmake_args[@]}" \
        -DCMAKE_MAKE_PROGRAM="$SYSTEM_NINJA" \
        -DCMAKE_C_COMPILER="$SYSTEM_CC" -DCMAKE_CXX_COMPILER="$SYSTEM_CXX" \
        -DCMAKE_BUILD_TYPE=Release -DBUILD_CONDA=OFF \
        -DSIPROS_BUILD_MPI=OFF \
        -DAERITH_ENABLE_TORCH=ON -DAERITH_TORCH_CPU_ONLY=ON \
        -DAERITH_TORCH_ROOT="$SYSTEM_TORCH_ROOT" \
        -DHDF5_USE_STATIC_LIBRARIES=OFF \
        -DHDF5_DIR="$RELEASE_HDF5_DIR" "$REPO_DIR"
    "$SYSTEM_NINJA" sipros aerith
    "$RELEASE_PREFIX/bin/patchelf" --force-rpath --set-rpath '$ORIGIN/lib' \
        "$SYSTEM_BUILD_DIR/bin/sipros"
    "$RELEASE_PREFIX/bin/patchelf" --force-rpath --set-rpath '$ORIGIN/lib' \
        "$SYSTEM_BUILD_DIR/bin/aerith"
    verify_release_dynamic_binary "$SYSTEM_BUILD_DIR/bin/sipros" \
        'libgomp\.so' 'libhdf5(_cpp)?\.so' 'libstdc\+\+\.so'
    verify_release_dynamic_binary "$SYSTEM_BUILD_DIR/bin/aerith" \
        'libtorch_cpu\.so' 'libc10\.so' 'libgomp\.so' \
        'libhdf5(_cpp)?\.so' 'libstdc\+\+\.so'
    
    # copy repo-built runtime commands atomically for publish/workflow use
    stage_publish_tools "$SYSTEM_BUILD_DIR/bin" cpu-release
    verify_glibc_217 "$REPO_DIR/tools/sipros"
    verify_packaged_release_binary "$REPO_DIR/tools/sipros" "$REPO_DIR/tools/lib"
    verify_packaged_release_binary "$REPO_DIR/tools/aerith" "$REPO_DIR/tools/lib" 1
    verify_glibc_217 "$REPO_DIR/tools/aerith" "$REPO_DIR/tools/lib"
    ;;
"buildConda")
    export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
    eval "$("$MAMBA_EXE" shell hook --shell=bash)"
    micromamba activate sipros5
    resolve_torch_cmake
    configure_mpi
    saved_ld_library_path="${LD_LIBRARY_PATH-}"
    unset LD_LIBRARY_PATH
    prepare_conda_build_dir "$CONDA_BUILD_DIR"
    conda_cuda_args=()
    if [ -x "$CONDA_PREFIX/bin/nvcc" ]; then
        conda_cuda_args+=(
            "-DCMAKE_CUDA_COMPILER=$CONDA_PREFIX/bin/nvcc"
            "-DCUDAToolkit_ROOT=$CONDA_PREFIX"
        )
    fi
    cd "$CONDA_BUILD_DIR"
    cmake -G Ninja "${cmake_args[@]}" -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_CONDA=ON \
        -DAERITH_ENABLE_TORCH=ON -DTorch_DIR="$TORCH_CMAKE_DIR" \
        -DAERITH_TORCH_PACKAGE_VERSION="$TORCH_VERSION" \
        "${conda_cuda_args[@]}" \
        -DHDF5_USE_STATIC_LIBRARIES=OFF \
        -DHDF5_DIR="$CONDA_PREFIX/cmake" "$REPO_DIR"
    export LD_LIBRARY_PATH="${saved_ld_library_path:+$saved_ld_library_path:}${CONDA_PREFIX}/lib"
    ninja
    verify_fully_dynamic_conda "$CONDA_BUILD_DIR/bin"
    verify_dynamic_torch "$CONDA_BUILD_DIR/bin/aerith"
    # add share lib for mpi version
    cd ..
    # deplist=$(ldd bin/siprosMPI | awk '{if (match($3,"/")){ print $3}}')
    # mkdir bin/libSiprosMPI
    # cp -L -n $deplist bin/libSiprosMPI

    # copy repo-built runtime commands atomically for publish/workflow use
    stage_publish_tools "$CONDA_BUILD_DIR/bin"
    ;;
"debug")
    export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
    eval "$($MAMBA_EXE shell hook --shell=bash)"
    micromamba activate sipros5
    resolve_torch_cmake
    configure_mpi
    saved_ld_library_path="${LD_LIBRARY_PATH-}"
    unset LD_LIBRARY_PATH
    mkdir -p "$CONDA_DEBUG_BUILD_DIR"
    cd "$CONDA_DEBUG_BUILD_DIR"
    cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS_DEBUG="-O0 -g3" \
        -DBUILD_CONDA=true -DAERITH_ENABLE_TORCH=ON \
        -DTorch_DIR="$TORCH_CMAKE_DIR" "$REPO_DIR"
    export LD_LIBRARY_PATH="${saved_ld_library_path:+$saved_ld_library_path:}${CONDA_PREFIX}/lib"
    ninja
    ;;
"make")
    "$SYSTEM_NINJA" -C "$SYSTEM_BUILD_DIR"
    ;;
"package")
    # Run clean and build before packaging
    $0 clean
    $0 build
    tmpdir=$(mktemp -d)
    trap 'rm -rf "$tmpdir"' EXIT
    mkdir -p "$tmpdir/sipros"
    cp -a tools script33 LICENSE "$tmpdir/sipros"
    rm -f "$tmpdir/sipros/tools/siprosMPI"
    rm -rf "$tmpdir/sipros/tools/lib"
    stage_release_runtime "$tmpdir/sipros/tools/sipros" \
        "$tmpdir/sipros/tools/lib"
    stage_release_runtime "$tmpdir/sipros/tools/aerith" \
        "$tmpdir/sipros/tools/lib"
    verify_glibc_217 "$tmpdir/sipros/tools/sipros"
    verify_packaged_release_binary "$tmpdir/sipros/tools/sipros" \
        "$tmpdir/sipros/tools/lib"
    verify_packaged_release_binary "$tmpdir/sipros/tools/aerith" \
        "$tmpdir/sipros/tools/lib" 1
    verify_glibc_217 "$tmpdir/sipros/tools/aerith" \
        "$tmpdir/sipros/tools/lib"
    if [ -f siprosRelease.zip ]; then
        rm siprosRelease.zip
    fi
    cd "$tmpdir"
    zip -r "$OLDPWD/siprosRelease.zip" "sipros" \
        -x "sipros/script33/debugProxy.py" \
        -x "*/__pycache__/*"
    cd "$OLDPWD"
    rm -rf "$tmpdir"
    trap - EXIT
    echo "Package created: siprosRelease.zip"
    ;;
"run")
    echo "Use the HDF5 workflow entrypoint, for example:"
    echo "  python script33/main.py -i data/pct1/raw/Pan_062822_X1iso5.raw -f data/EcoliWithCrapNodup.fasta -o data/tmp/raxport_hdf5_workflow_test/direct_fasta -t 4"
    ;;
*)
    ./make "build"
    ;;
esac
