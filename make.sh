#!/bin/bash
set -e
## install compile dependencies with apt
# sudo apt install python openmpi-bin libopenmpi-dev build-essential cmake ninja-build gdb google-perftools libgoogle-perftools-dev
## install compile dependencies from micromamba 
# micromamba create -n sipros5 -c conda-forge openmpi gxx_linux-64 gcc_linux-64 cmake ninja gdb gperftools python=3.12 lxml pandas
## compiler name x86_64-conda_cos6-linux-gnu-g++
# micromamba activate sipros5
## run follows to load dynamic libs when running bin/SiprosV3omp bin/SiprosV3mpi bin/SiprosV3test
# micromamba activate sipros5
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAMBA_EXE="${MAMBA_EXE:-micromamba}"
if ! command -v "$MAMBA_EXE" >/dev/null 2>&1 && [ -x "$HOME/.local/bin/micromamba" ]; then
    MAMBA_EXE="$HOME/.local/bin/micromamba"
fi
VCPKG_ROOT="${VCPKG_ROOT:-$REPO_DIR/vcpkg}"
VCPKG_TARGET_TRIPLET="${VCPKG_TARGET_TRIPLET:-x64-linux}"
VCPKG_TOOLCHAIN_FILE="${VCPKG_TOOLCHAIN_FILE:-$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake}"

require_vcpkg_toolchain() {
    if [ ! -f "$VCPKG_TOOLCHAIN_FILE" ]; then
        echo "Missing vcpkg toolchain: $VCPKG_TOOLCHAIN_FILE" >&2
        echo "Install vcpkg at $VCPKG_ROOT first." >&2
        exit 1
    fi
}

prepare_vcpkg_build_dir() {
    mkdir -p "$1"
    local cache="$1/CMakeCache.txt"
    local hdf5_dir="HDF5_DIR:PATH=$VCPKG_ROOT/installed/$VCPKG_TARGET_TRIPLET/share/hdf5"
    if [ -f "$cache" ] && ! grep -Fq "$hdf5_dir" "$cache"; then
        echo "Reconfiguring $1 with vcpkg toolchain"
        rm -rf "$cache" "$1/CMakeFiles"
    fi
}

prepare_conda_build_dir() {
    mkdir -p "$1"
    local cache="$1/CMakeCache.txt"
    local hdf5_dir="HDF5_DIR:PATH=$CONDA_PREFIX/cmake"
    if [ -f "$cache" ] && ! grep -Fq "$hdf5_dir" "$cache"; then
        echo "Reconfiguring $1 with micromamba HDF5"
        rm -rf "$cache" "$1/CMakeFiles"
    fi
}

run_sipros5() {
    "$MAMBA_EXE" run -n sipros5 "$@"
}

cmake_args=()
if [ -n "${CMAKE_ARGS:-}" ]; then
    cmake_args=(${CMAKE_ARGS})
fi

case $1 in
"load") ;;
"clean")
    rm -rf build
    mkdir build
    rm -rf bin
    mkdir bin
    ;;
"build")
    require_vcpkg_toolchain
    prepare_vcpkg_build_dir build
    mkdir -p tools
    cd build
    run_sipros5 cmake -G Ninja "${cmake_args[@]}" -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE="$VCPKG_TOOLCHAIN_FILE" -DVCPKG_TARGET_TRIPLET="$VCPKG_TARGET_TRIPLET" ..
    run_sipros5 ninja
    # add share lib for mpi version
    cd ..
    # deplist=$(ldd bin/siprosMPI | awk '{if (match($3,"/")){ print $3}}')
    # mkdir bin/libSiprosMPI
    # cp -L -n $deplist bin/libSiprosMPI
    
    # copy file atomic when running
    tmpdir=$(mktemp -d tools/.tmp.XXXXXX)
    cp bin/* "$tmpdir"
    mv "$tmpdir"/* tools/
    rmdir "$tmpdir"
    ;;
"buildConda")
    export MAMBA_ROOT_PREFIX=~/micromamba
    eval "$("$MAMBA_EXE" shell hook --shell=bash)"
    micromamba activate sipros5
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
    prepare_conda_build_dir build
    cd build
    cmake -G Ninja "${cmake_args[@]}" -DCMAKE_BUILD_TYPE=Release -DBUILD_CONDA=true -DHDF5_DIR="$CONDA_PREFIX/cmake" ..
    ninja
    # add share lib for mpi version
    cd ..
    # deplist=$(ldd bin/siprosMPI | awk '{if (match($3,"/")){ print $3}}')
    # mkdir bin/libSiprosMPI
    # cp -L -n $deplist bin/libSiprosMPI

    # copy file atomic when running
    tmpdir=$(mktemp -d tools/.tmp.XXXXXX)
    cp bin/* "$tmpdir"
    mv "$tmpdir"/* tools/
    rmdir "$tmpdir"
    ;;
"buildTick")
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DTicktock=Ticktock ..
    make -j8
    ;;
"debug")
    export MAMBA_ROOT_PREFIX=~/micromamba
    eval "$(~/.local/bin/micromamba shell hook --shell=bash)"
    micromamba activate sipros5
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
    cd build
    cmake -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS_DEBUG="-O0 -g3" -DBUILD_CONDA=true ..
    ninja
    ;;
"make")
    cd build
    make .. -j8
    ;;
"package")
    # Run clean and build before packaging
    $0 clean
    $0 build
    tmpdir=$(mktemp -d)
    mkdir -p "$tmpdir/sipros"
    cp -r configTemplates tools script33 LICENSE "$tmpdir/sipros"
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
    cd timeCompare
    starttime=$(date +'%Y-%m-%d %H:%M:%S')
    ../bin/SiprosV3omp -c ../../data/SiproConfig.N15_0Pct.cfg \
        -f ../../data/AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 \
        -o . -s
    endtime=$(date +'%Y-%m-%d %H:%M:%S')
    start_seconds=$(date --date="$starttime" +%s)
    end_seconds=$(date --date="$endtime" +%s)
    echo "running time： "$((end_seconds - start_seconds))"s"
    ;;
*)
    ./make "build"
    ;;
esac
