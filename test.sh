#!/bin/bash
source ~/miniconda3/etc/profile.d/conda.sh
binPath="/ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP"
cd test
case $1 in
"mkdir")
    mkdir raw ft mzml sip fasta regular configs
    /scratch/yixiong/benchmark/configGenerator -i SiprosEnsembleC13.cfg -o configs -e C
    ;;
"mkdb")
    # conda activate py2
    conda activate sklearn
    # python ../script/sipros_prepare_protein_database.py \
    python ../script3/sipros_prepare_protein_database.py \
        -i fasta/EcoliWithCrapNodup.fasta \
        -o fasta/Decoy.fasta \
        -c SiprosEnsembleConfig.cfg
    python ../script3/sipros_prepare_protein_database.py \
        -i fasta/EcoliWithCrapNodup.fasta \
        -o fasta/DecoySIP.fasta \
        -c SiprosEnsembleC13.cfg
    ${binPath}/configGenerator \
        -i SiprosEnsembleC13.cfg \
        -o configs \
        -e C
    ;;
"convert")
    conda activate mono
    mono ${binPath}/Raxport.exe -i raw -o ft &
    # ensemble only support indexed mzml
    mono tools/ThermoRawFileParser/ThermoRawFileParser.exe -d=raw -f=1 -o mzml
    ;;
"clean")
    rm -r regular/* sip/*
    ;;
"run")
    printf "\n=====Regular Search=====\n\n"
    conda activate mpi
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
    starttime=$(date +'%Y-%m-%d %H:%M:%S')
    # ${binPath}/SiprosEnsembleOMP -w ft -c SiprosEnsembleConfig.cfg -o regular
    # ../bin/SiprosEnsembleOMP -w mzml -c SiprosEnsembleConfig.cfg -o regular
    # ../bin/SiprosEnsembleOMP -w ft -c SiprosEnsembleConfig.cfg -o regular
    # /ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP/SiprosEnsembleOMP -w mzml -c SiprosEnsembleConfig.cfg -o regular
    # ft='/ourdisk/hpc/prebiotics/yixiong/auto_archive_notyet/ubuntuShare/EcoliSIP/goodResults/pct1ensemble/ft'
    ft='/ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/UbuntuShare/benchmark/pct1/ft'
    ../bin/SiprosEnsembleOMP -w ${ft} -c SiprosEnsembleConfig.cfg -o regular
    endtime=$(date +'%Y-%m-%d %H:%M:%S')
    start_seconds=$(date --date="$starttime" +%s)
    end_seconds=$(date --date="$endtime" +%s)
    echo "running time： "$((end_seconds - start_seconds))"s"
    wait

    printf "\n=====Filter PSM=====\n\n"
    conda activate sklearn
    starttime=$(date +'%Y-%m-%d %H:%M:%S')
    python ../script3/sipros_psm_tabulating.py \
        -i regular \
        -c SiprosEnsembleConfig.cfg -o regular
    python ../script3/sipros_ensemble_filtering.py \
        -i regular \
        -c SiprosEnsembleConfig.cfg \
        -o regular
    python ../script3/sipros_peptides_assembling.py \
        -c SiprosEnsembleConfig.cfg \
        -w regular
    endtime=$(date +'%Y-%m-%d %H:%M:%S')
    start_seconds=$(date --date="$starttime" +%s)
    end_seconds=$(date --date="$endtime" +%s)
    echo "running time： "$((end_seconds - start_seconds))"s"
    ;;
"runSIPone")
    conda activate mpi
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
    ../bin/SiprosEnsembleOMP -w ft -c configs/C13_5000Pct.cfg -o sip
    ;;
"runSIP")
    printf "\n=====SIP Search=====\n\n"
    conda activate mpi
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${CONDA_PREFIX}/lib
    starttime=$(date +'%Y-%m-%d %H:%M:%S')

    export OMP_NUM_THREADS=10
    export configs=(configs/*.cfg)
    export ft='/ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/UbuntuShare/benchmark/pct1/ft'
    # export ft='/ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/UbuntuShare/benchmark/pct50/ft'
    # export ft='/ourdisk/hpc/nullspace/yixiong/auto_archive_notyet/tape_2copies/UbuntuShare/benchmark/pct99/ft'
    echo "${configs[@]}" | xargs -n 1 -P 10 \
        bash -c '../bin/SiprosEnsembleOMP -w ${ft} -c $0 -o sip'

    endtime=$(date +'%Y-%m-%d %H:%M:%S')
    start_seconds=$(date --date="$starttime" +%s)
    end_seconds=$(date --date="$endtime" +%s)
    echo "running time： "$((end_seconds - start_seconds))"s"
    wait

    printf "\n=====Filter PSM=====\n\n"
    conda activate sklearn
    starttime=$(date +'%Y-%m-%d %H:%M:%S')
    python ../script3/sipros_psm_tabulating.py \
        -i sip \
        -c SiprosEnsembleC13.cfg -o sip
    python ../script3/sipros_ensemble_filtering.py \
        -i sip \
        -c SiprosEnsembleC13.cfg \
        -o sip
    python ../script3/sipros_peptides_assembling.py \
        -c SiprosEnsembleC13.cfg \
        -w sip
    endtime=$(date +'%Y-%m-%d %H:%M:%S')
    start_seconds=$(date --date="$starttime" +%s)
    end_seconds=$(date --date="$endtime" +%s)
    echo "running time： "$((end_seconds - start_seconds))"s"
    ;;
"filter")
    conda activate pypy2
    printf "\n=====Filter PSM=====\n\n"
    python ../script/sipros_psm_tabulating.py \
        -i regular \
        -c SiprosEnsembleConfig.cfg -o regular
    python ../script/sipros_ensemble_filtering.py \
        -i regular/*.tab \
        -c SiprosEnsembleConfig.cfg \
        -o regular
    printf "\n====Assemble protein=====\n\n"
    python ../script/sipros_peptides_assembling.py \
        -c SiprosEnsembleConfig.cfg \
        -w regular
    ;;
*)
    ./cmd.sh "run"
    ;;
esac
