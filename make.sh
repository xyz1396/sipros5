#!/bin/bash
set -e

# Build dependency setup
#
# System tools used by ./make.sh build:
#   sudo apt update
#   sudo apt install build-essential git curl cmake ninja-build openmpi-bin libopenmpi-dev
#
# set sysroot_linux-64=2.17 for centos7, 2.34 for rocky9
# 
# Micromamba/Conda build (./make.sh buildConda):
#   micromamba create -n sipros5 -c conda-forge \
#     hdf5 h5py openmpi cmake ninja gcc_linux-64 gxx_linux-64 gdb gperftools \
#     python=3.12 lxml pandas sysroot_linux-64=2.34 matplotlib
#
# Static HDF5 used by the system build:
#   git clone --depth 1 https://github.com/microsoft/vcpkg.git vcpkg
#   ./vcpkg/bootstrap-vcpkg.sh -disableMetrics
#   ./vcpkg/vcpkg install "hdf5[cpp]:x64-linux"
#
# Activate the environment before running Conda-built binaries so their
# OpenMP, compiler-runtime, and MPI shared libraries are available:
#   micromamba activate sipros5
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAMBA_EXE="${MAMBA_EXE:-micromamba}"
if ! command -v "$MAMBA_EXE" >/dev/null 2>&1 && [ -x "$HOME/.local/bin/micromamba" ]; then
    MAMBA_EXE="$HOME/.local/bin/micromamba"
fi
VCPKG_ROOT="${VCPKG_ROOT:-$REPO_DIR/vcpkg}"
VCPKG_TARGET_TRIPLET="${VCPKG_TARGET_TRIPLET:-x64-linux}"
VCPKG_TOOLCHAIN_FILE="${VCPKG_TOOLCHAIN_FILE:-$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake}"
VCPKG_HDF5_DIR="${VCPKG_HDF5_DIR:-$VCPKG_ROOT/installed/$VCPKG_TARGET_TRIPLET/share/hdf5}"
SYSTEM_CMAKE="${SYSTEM_CMAKE:-/usr/bin/cmake}"
SYSTEM_NINJA="${SYSTEM_NINJA:-/usr/bin/ninja}"
SYSTEM_CC="${SYSTEM_CC:-/usr/bin/gcc}"
SYSTEM_CXX="${SYSTEM_CXX:-/usr/bin/g++}"
SYSTEM_MPI_CXX="${SYSTEM_MPI_CXX:-/usr/bin/mpicxx}"
BUILD_ROOT="$REPO_DIR/build"
SYSTEM_BUILD_DIR="$BUILD_ROOT/system"
CONDA_BUILD_DIR="$BUILD_ROOT/conda"
CONDA_DEBUG_BUILD_DIR="$BUILD_ROOT/conda-debug"
RUNTIME_TOOL_BINARIES=(
    Raxport-linux-x64
    philosopher-v5.1.2
    sipros
    siprosMPI
    aerith
    deepfilter
    metaLP
    ionquant
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

require_vcpkg_toolchain() {
    if [ ! -f "$VCPKG_TOOLCHAIN_FILE" ]; then
        echo "Missing vcpkg toolchain: $VCPKG_TOOLCHAIN_FILE" >&2
        echo "Install vcpkg at $VCPKG_ROOT first." >&2
        exit 1
    fi
}

require_executable() {
    if [ ! -x "$1" ]; then
        echo "Missing required executable: $1" >&2
        exit 1
    fi
}

require_static_hdf5() {
    local hdf5_dir="$1"
    local triplet="$2"
    if [ ! -f "$hdf5_dir/hdf5-config.cmake" ]; then
        echo "Missing vcpkg HDF5 CMake package: $hdf5_dir" >&2
        echo "Install hdf5[cpp]:$triplet with vcpkg first." >&2
        exit 1
    fi
}

prepare_vcpkg_build_dir() {
    mkdir -p "$1"
    local cache="$1/CMakeCache.txt"
    local hdf5_dir="=$VCPKG_HDF5_DIR"
    local compiler="=$SYSTEM_CXX"
    if [ -f "$cache" ] &&
       { ! grep -F 'HDF5_DIR:' "$cache" | grep -Fq "$hdf5_dir" ||
         ! grep -F 'CMAKE_CXX_COMPILER:' "$cache" | grep -Fq "$compiler"; }; then
        echo "Reconfiguring $1 with vcpkg toolchain"
        rm -rf "$cache" "$1/CMakeFiles"
    fi
}

prepare_conda_build_dir() {
    mkdir -p "$1"
    local cache="$1/CMakeCache.txt"
    local hdf5_dir="=$CONDA_PREFIX/cmake"
    if [ -f "$cache" ] &&
       ! grep -F 'HDF5_DIR:' "$cache" | grep -Fq "$hdf5_dir"; then
        echo "Reconfiguring $1 with dynamic Conda HDF5"
        rm -rf "$cache" "$1/CMakeFiles"
    fi
}

run_sipros5() {
    "$MAMBA_EXE" run -n sipros5 "$@"
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

verify_fully_static_omp() {
    local binary="$1/sipros"
    local dependencies
    dependencies=$(ldd "$binary" 2>&1 || true)
    if [[ "$dependencies" != *"not a dynamic executable"* &&
          "$dependencies" != *"statically linked"* ]]; then
        echo "Expected a fully static OpenMP binary: $binary" >&2
        echo "$dependencies" >&2
        return 1
    fi
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

stage_publish_tools() {
    local source_bin_dir="$1"
    local binaries=(sipros siprosMPI aerith)
    local destinations=("$REPO_DIR/bin" "$REPO_DIR/tools")
    local binary destination tmpdir

    # Validate the complete build before replacing any published executable.
    for binary in "${binaries[@]}"; do
        if [ ! -x "$source_bin_dir/$binary" ]; then
            echo "Missing built executable: $source_bin_dir/$binary" >&2
            return 1
        fi
    done

    for destination in "${destinations[@]}"; do
        mkdir -p "$destination"
        tmpdir=$(mktemp -d "$destination/.stage.XXXXXX")
        for binary in "${binaries[@]}"; do
            install -m 0755 "$source_bin_dir/$binary" "$tmpdir/$binary"
        done
        for binary in "${binaries[@]}"; do
            mv -f "$tmpdir/$binary" "$destination/$binary"
        done
        rmdir "$tmpdir"
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
    ;;
"build")
    require_vcpkg_toolchain
    require_static_hdf5 "$VCPKG_HDF5_DIR" "$VCPKG_TARGET_TRIPLET"
    require_executable "$SYSTEM_CMAKE"
    require_executable "$SYSTEM_NINJA"
    require_executable "$SYSTEM_CC"
    require_executable "$SYSTEM_CXX"
    require_executable "$SYSTEM_MPI_CXX"
    cmake_args+=("-DMPI_CXX_COMPILER=$SYSTEM_MPI_CXX" "-DMPI_CXX_SKIP_MPICXX=ON")
    prepare_vcpkg_build_dir "$SYSTEM_BUILD_DIR"
    mkdir -p tools
    cd "$SYSTEM_BUILD_DIR"
    CC="$SYSTEM_CC" CXX="$SYSTEM_CXX" "$SYSTEM_CMAKE" -G Ninja "${cmake_args[@]}" \
        -DCMAKE_MAKE_PROGRAM="$SYSTEM_NINJA" \
        -DCMAKE_C_COMPILER="$SYSTEM_CC" -DCMAKE_CXX_COMPILER="$SYSTEM_CXX" \
        -DCMAKE_BUILD_TYPE=Release -DSIPROS_STATIC_HDF5=ON \
        -DCMAKE_TOOLCHAIN_FILE="$VCPKG_TOOLCHAIN_FILE" \
        -DVCPKG_TARGET_TRIPLET="$VCPKG_TARGET_TRIPLET" \
        -DZLIB_DIR="$VCPKG_ROOT/installed/$VCPKG_TARGET_TRIPLET/share/zlib" \
        -Dlibaec_DIR="$VCPKG_ROOT/installed/$VCPKG_TARGET_TRIPLET/share/libaec" \
        -DHDF5_DIR="$VCPKG_HDF5_DIR" "$REPO_DIR"
    "$SYSTEM_NINJA"
    verify_fully_static_omp "$SYSTEM_BUILD_DIR/bin"
    # add share lib for mpi version
    cd ..
    # deplist=$(ldd bin/siprosMPI | awk '{if (match($3,"/")){ print $3}}')
    # mkdir bin/libSiprosMPI
    # cp -L -n $deplist bin/libSiprosMPI
    
    # copy repo-built runtime commands atomically for publish/workflow use
    stage_publish_tools "$SYSTEM_BUILD_DIR/bin"
    ;;
"buildConda")
    export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
    eval "$("$MAMBA_EXE" shell hook --shell=bash)"
    micromamba activate sipros5
    configure_mpi
    saved_ld_library_path="${LD_LIBRARY_PATH-}"
    unset LD_LIBRARY_PATH
    prepare_conda_build_dir "$CONDA_BUILD_DIR"
    cd "$CONDA_BUILD_DIR"
    cmake -G Ninja "${cmake_args[@]}" -DCMAKE_BUILD_TYPE=Release \
        -DBUILD_CONDA=ON -DSIPROS_STATIC_HDF5=OFF \
        -DHDF5_USE_STATIC_LIBRARIES=OFF \
        -DHDF5_DIR="$CONDA_PREFIX/cmake" "$REPO_DIR"
    export LD_LIBRARY_PATH="${saved_ld_library_path:+$saved_ld_library_path:}${CONDA_PREFIX}/lib"
    ninja
    verify_fully_dynamic_conda "$CONDA_BUILD_DIR/bin"
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
    configure_mpi
    saved_ld_library_path="${LD_LIBRARY_PATH-}"
    unset LD_LIBRARY_PATH
    mkdir -p "$CONDA_DEBUG_BUILD_DIR"
    cd "$CONDA_DEBUG_BUILD_DIR"
    cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS_DEBUG="-O0 -g3" \
        -DBUILD_CONDA=true "$REPO_DIR"
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
    mkdir -p "$tmpdir/sipros"
    cp -a tools script33 LICENSE "$tmpdir/sipros"
    if [ -f siprosRelease.zip ]; then
        rm siprosRelease.zip
    fi
    cd "$tmpdir"
    zip -r "$OLDPWD/siprosRelease.zip" "sipros" \
        -x "sipros/script33/debugProxy.py" "sipros/script33/quant.py" \
        -x "*/__pycache__/*"
    cd "$OLDPWD"
    rm -rf "$tmpdir"
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
