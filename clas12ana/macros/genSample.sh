#!/bin/bash
###########################################################################################
###########################################################################################
###########################################################################################
PROJECT_NAME="pipluspi0_noresonance"
NNPI0DIR=/work/clas12/users/gmat/nnPi0
###########################################################################################
###########################################################################################
###########################################################################################
USERNAME="$USER"
pwd=$PWD
###########################################################################################
###########################################################################################
###########################################################################################
printred () {
  echo -e "\e[31m$1\e[0m"
}

printblue () {
  echo -e "\e[34m$1\e[0m"
}

printgreen () {
  echo -e "\e[32m$1\e[0m"
}
###########################################################################################
###########################################################################################
###########################################################################################
hl="--------------------------------------------------------------------------------------------------------------------------------------------"

if [ $# -lt 1 ]; then
  echo """
  USAGE: $0 [ML Method] [flags(optional)]
  Automates the generation of raw TTrees from the CLAS12 HIPO data
  After the raw TTrees are created (or if they are created already) this program preprocesses
  the raw TTrees into the input space of various ML algorithms
   - [ML Method]: Options (catboost particleNet)            
   - [flags]:
              -o      (overwrite existing raw TTrees)
              -nFiles <INT> (maximum number of files per hipo directory: default=5)
              -noslurm (run all the processing code in command line without slurm)
              -maxEvents <INT> (set the maximum number of events analyzed per hipo file)
  """
  exit 2
fi

MLmethods="catboost particleNet"
nCPUs=4
memPerCPU=4000

declare -A flags
declare -A booleans
args=()

while [ "$1" ];
do
    arg=$1
    if [ "${1:0:1}" == "-" ]
    then
      shift
      rev=$(echo "$arg" | rev)
      if [ -z "$1" ] || [ "${1:0:1}" == "-" ] || [ "${rev:0:1}" == ":" ]
      then
        bool=$(echo ${arg:1} | sed s/://g)
        booleans[$bool]=true
      else
        value=$1
        flags[${arg:1}]=$value
        shift
      fi
    else
      args+=("$arg")
      shift
    fi
done

#######################################################################################################
MLmethod=${args[0]}
nFiles=5

if [ ! -z ${flags["nFiles"]} ]; then
    nFiles=${flags["nFiles"]}
fi

if ! echo "$MLmethods" | grep -q "$MLmethod"; then
    echo $hl
    echo "ERROR: ML Method $MLmethod not valid. Must use from following list [$MLmethods]"
    exit 2
else
    rungroup="${rungroup//-}"
fi

nEvents=1000000 # per file
if [ ! -z ${flags["maxEvents"]} ]; then
    nEvents=${flags["maxEvents"]}
fi

#######################################################################################################
#######################################################################################################
#######################################################################################################
# CREATE PROJECT DIRECTORY 
# (this is created in BOTH /work and /volatile)
#######################################################################################################
#######################################################################################################
#######################################################################################################
rawdir=/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/raw
volatiledir=/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/projects/$PROJECT_NAME
workdir=/work/clas12/users/gmat/scipio/projects/$PROJECT_NAME
printred "************ Creating /volatile directories for PROJECT_NAME=$PROJECT_NAME ************"
if [ -d "${volatiledir}/" ]; then
    printred "\t ${volatiledir} already exists...continuing..."
else
    printgreen "\t mkdir ${volatiledir}"
    mkdir -p ${volatiledir}
fi
if [ -d "${volatiledir}/${MLmethod}" ]; then
    printred "\t ${volatiledir}/${MLmethod} already exists...rewriting..."
    rm -r ${volatiledir}/${MLmethod}
fi
printgreen "\t mkdir ${volatiledir}/${MLmethod}"
mkdir -p ${volatiledir}/${MLmethod}
printgreen "\t mkdir ${volatiledir}/${MLmethod}/MLinput"
mkdir -p ${volatiledir}/${MLmethod}/MLinput
printgreen "\t mkdir ${volatiledir}/${MLmethod}/MLoutput"
mkdir -p ${volatiledir}/${MLmethod}/MLoutput
printgreen "\t mkdir ${volatiledir}/${MLmethod}/postprocess"
mkdir -p ${volatiledir}/${MLmethod}/postprocess
printgreen "\t mkdir ${volatiledir}/${MLmethod}/postprocess_binned"
mkdir -p ${volatiledir}/${MLmethod}/postprocess_binned

printred "************Creating /work directories for PROJECT_NAME=$PROJECT_NAME************"
if [ -d "${workdir}/" ]; then
    printred "\t ${workdir} already exists...continuing..."
else
    printgreen "\t mkdir ${workdir}"
    mkdir -p $workdir
fi
if [ -d "${workdir}/${MLmethod}" ]; then
    printred "\t ${workdir}/${MLmethod} already exists...rewriting..."
    rm -r ${workdir}/${MLmethod}
fi
printgreen "\t mkdir ${workdir}/${MLmethod}"
mkdir -p $workdir/$MLmethod

#######################################################################################################
#######################################################################################################
#######################################################################################################
# In the ${workdir}/${MLmethod} subdirectory, create a training and predicting script for the ML
if [ $MLmethod == "catboost" ]; then
    # For catboost, we need the input spreadsheet, a train.sh script, and a predict.sh script
    # For now, just copy the following default input yaml
    CATBOOST_INPUT_YAML=$workdir/$MLmethod/MLinputs.yaml
    CATBOOST_MODEL_INBENDING=$workdir/$MLmethod/model_inbending
    CATBOOST_MODEL_OUTBENDING=$workdir/$MLmethod/model_outbending
    printgreen "\t Generating MLinputs.yaml file in $workdir/$MLmethod"
    cp $NNPI0DIR/catboost/input/input_noresonance.yaml $workdir/$MLmethod/MLinputs.yaml
    
    # Create the train.sh script
    printgreen "\t Generating train.sh file in $workdir/$MLmethod"
    cat >> $workdir/$MLmethod/train.sh << EOF
cd $NNPI0DIR/catboost/
python train.py --data_dir ${volatiledir}/${MLmethod}/MLinput \
--depth 12 \
--make_plots true \
--subdata MC_inbending \
--model_dir $CATBOOST_MODEL_INBENDING \
--input_yaml $CATBOOST_INPUT_YAML

python train.py --data_dir ${volatiledir}/${MLmethod}/MLinput \
--depth 12 \
--make_plots true \
--subdata MC_outbending \
--model_dir $CATBOOST_MODEL_OUTBENDING \
--input_yaml $CATBOOST_INPUT_YAML

cd $pwd
EOF
    # Create the predict.sh script
    printgreen "\t Generating predict.sh file in $workdir/$MLmethod"
    cat >> $workdir/$MLmethod/predict.sh << EOF
cd $NNPI0DIR/catboost/
python predict.py \
    --version nSidis \
    --model_dir $CATBOOST_MODEL_INBENDING \
    --subdata RGA_inbending \
    --output_dir ${volatiledir}/${MLmethod}/MLoutput
    
python predict.py \
    --version nSidis \
    --model_dir $CATBOOST_MODEL_OUTBENDING \
    --subdata RGA_outbending \
    --output_dir ${volatiledir}/${MLmethod}/MLoutput
    
python predict.py \
    --version MC \
    --model_dir $CATBOOST_MODEL_INBENDING \
    --subdata MC_inbending \
    --output_dir ${volatiledir}/${MLmethod}/MLoutput
    
python predict.py \
    --version MC \
    --model_dir $CATBOOST_MODEL_OUTBENDING \
    --subdata MC_outbending \
    --output_dir ${volatiledir}/${MLmethod}/MLoutput
    
cd $pwd
EOF


fi

#######################################################################################################
#######################################################################################################
#######################################################################################################


slurm=true
if [ ! -z ${booleans["noslurm"]} ]; then
    slurm=false
fi



logdir=""
slurmdir=""
if [ $slurm == true ]; then 
    
        NOW=$( date '+%F_%H_%M_%S' )
        NOWdir=/farm_out/gmat/clas12analysis.sidis.data/rga/ML/$NOW
        if [ -d "${NOWdir}/" ]; then
            rm -r ${NOWdir}
        fi
        mkdir ${NOWdir}
        logdir=$NOWdir"/log"
        slurmdir=$NOWdir"/slurm"
        if [ -d "${logdir}/" ]; then
            rm -r ${logdir}
        fi
        if [ -d "${slurmdir}/" ]; then
            rm -r ${slurmdir}
        fi
        mkdir ${logdir}
        mkdir ${slurmdir}
fi



nCPUs=4
memPerCPU=4000
versions=("nSidis" "MC")
for ana in "${versions[@]}"
do
    printblue "\n\nVERSION=$ana\n\n"
    hipo_is_mc=0
    if [ $ana == "MC" ]; then
        declare -a hipofiles=("/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg*/*.hipo" "/cache/clas12/rg-a/production/montecarlo/clasdis/fall\2018/torus+1/v1/bkg*/*.hipo")
        hipo_is_mc=1
    elif [ $ana == "nSidis" ]; then
        declare -a hipofiles=("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis/*.hipo" "/cache/clas12/rg-a/production/recon/fall2018/tor\us+1/pass1/v1/dst/train/nSidis/*.hipo" "/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass1/v1/dst/train/nSidis/*.hipo")
    fi



    for subhipofiles in "${hipofiles[@]}"
    do
    filenum=0
    for hipo in $subhipofiles
    do
        beamE=0
        runNumber=99999
        if [ $ana == "MC" ]; then
        beamE=10.604
        read runNumber <<< $(basename $hipo | grep -oP '(?<=_job_).*(?=.hipo)')
        else
        read runNumber <<< $(basename $hipo | grep -oP '(?<='$ana'_00).*(?=.hipo)')
            if echo $hipo | grep -w -q "spring2019"; then
                beamE=10.2
            else
                beamE=10.6
            fi
        fi

        raw_out=$rawdir/${ana}_${runNumber}.root
        MLmethod_out=$volatiledir"/$MLmethod/MLinput/${ana}_${runNumber}.root"

        if [ $slurm == true ]; then 

            slurmshell=${slurmdir}"/${ana}_${runNumber}_${MLmethod}.sh"
            slurmslurm=${slurmdir}"/${ana}_${runNumber}_${MLmethod}.slurm"

            touch $slurmshell
            touch $slurmslurm
            chmod +x $slurmshell

            cat >> $slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=${memPerCPU}
#SBATCH --job-name=job_${ana}_${runNumber}_${MLmethod}
#SBATCH --cpus-per-task=${nCPUs}
#SBATCH --time=24:00:00
#SBATCH --output=${logdir}/${ana}_${runNumber}_${MLmethod}.out
#SBATCH --error=${logdir}/${ana}_${runNumber}_${MLmethod}.err
$slurmshell
EOF

            if [ ! -z ${booleans["o"]} ]; then
                cat >> $slurmshell << EOF
#!/bin/tcsh
source /group/clas12/packages/setup.csh
module load clas12/pro
clas12root -b -q $PWD/pi0_readHipo.C\(\"${hipo}\",\"${raw_out}\",$beamE,$nEvents,$hipo_is_mc\)
clas12root -b -q $PWD/pi0_preprocess_${MLmethod}.C\(\"${raw_out}\",\"${MLmethod_out}\"\)
echo "Done"
EOF
            else
                cat >> $slurmshell << EOF
#!/bin/tcsh
source /group/clas12/packages/setup.csh
module load clas12/pro
if ( ! -e "${raw_out}" ) then
    clas12root -b -q $PWD/pi0_readHipo.C\(\"${hipo}\",\"${raw_out}\",$beamE,$nEvents,$hipo_is_mc\)
endif
clas12root -b -q $PWD/pi0_preprocess_${MLmethod}.C\(\"${raw_out}\",\"${MLmethod_out}\"\)
echo "Done"
EOF
            fi

            sbatch $slurmslurm

        else    

            printblue $hl
            printblue "Reading hipo file $hipo"
            printblue $hl

        if [ -f "$raw_out" ]; then
                if [ ! -z ${booleans["o"]} ]; then	
            clas12root -b -q "pi0_readHipo.C(\"${hipo}\",\"${raw_out}\",$beamE,$nEvents,$hipo_is_mc)"
            else		
            printred "$raw_out exists...skipping..."
                fi
        else
            clas12root -b -q "pi0_readHipo.C(\"${hipo}\",\"${raw_out}\",$beamE,$nEvents,$hipo_is_mc)"
            fi
            printblue $hl
            printblue "Preprocessing hipo file $hipo"
            printblue $hl
            clas12root -b -q "pi0_preprocess_${MLmethod}.C(\"${raw_out}\",\"${MLmethod_out}\")"
            printblue $hl

        fi

        filenum=$((filenum+1))
        if [ $filenum -eq $nFiles ]; then
            break
        fi

    done
    done


    echo $hl
    printgreen "Location of log files: $logdir"
    printgreen "Location of slurm files: $slurmdir"
done