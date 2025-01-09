#!/bin/bash
# micromamba create -n sipros5 -c conda-forge openmpi gxx_linux-64 gcc_linux-64 cmake gperftools
# compiler name x86_64-conda_cos6-linux-gnu-g++
# micromamba activate sipros5
# run follows to load dynamic libs when running bin/SiprosV3omp bin/SiprosV3mpi bin/SiprosV3test
# micromamba activate sipros5
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
case $1 in
"load") ;;
"clean")
    rm -rf build
    mkdir build
    rm -rf bin
    mkdir bin
    ;;
"build")
    mkdir build
    cd build
    cmake -G Ninja -DCMAKE_BUILD_TYPE=Release ..
    ninja
    # add share lib for mpi version
    cd ..
    # deplist=$(ldd bin/siprosMPI | awk '{if (match($3,"/")){ print $3}}')
    # mkdir bin/libSiprosMPI
    # cp -L -n $deplist bin/libSiprosMPI
    cp bin/* tools
    ;;
"buildConda")
    export MAMBA_ROOT_PREFIX=~/micromamba
    eval "$(~/.local/bin/micromamba shell hook --shell=bash)"
    micromamba activate sipros5
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
    mkdir build
    cd build
    cmake -G Ninja ${CMAKE_ARGS} -DCMAKE_BUILD_TYPE=Release -DBUILD_CONDA=true ..
    ninja
    # add share lib for mpi version
    cd ..
    # deplist=$(ldd bin/siprosMPI | awk '{if (match($3,"/")){ print $3}}')
    # mkdir bin/libSiprosMPI
    # cp -L -n $deplist bin/libSiprosMPI
    cp bin/* tools
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
