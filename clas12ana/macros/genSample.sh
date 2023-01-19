#!/bin/bash
###########################################################################################
###########################################################################################
###########################################################################################
PROJECT_NAME="pipluspi0_prelim"
NNPI0DIR=/work/clas12/users/gmat/nnPi0
SCIPIODIR=/work/clas12/users/gmat/scipio
###########################################################################################
###########################################################################################
###########################################################################################
USERNAME="$USER"
pwd=$PWD
nCPUs=4
memPerCPU=4000
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

nEvents=1000000000 # per file
if [ ! -z ${flags["maxEvents"]} ]; then
    nEvents=${flags["maxEvents"]}
fi

slurm=true
if [ ! -z ${booleans["noslurm"]} ]; then
    slurm=false
fi

#######################################################################################################
#######################################################################################################
#######################################################################################################

logdir=""
slurmdir=""
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
workdir=$SCIPIODIR/projects/$PROJECT_NAME
printred "************ Creating /volatile directories for PROJECT_NAME=$PROJECT_NAME ************"
if [ -d "${volatiledir}/" ]; then
    printred "\t ${volatiledir} already exists...aborting..."
    #exit
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
    rm -rf ${workdir}/${MLmethod}
fi
printgreen "\t mkdir ${workdir}/${MLmethod}"
mkdir -p $workdir/$MLmethod
printgreen "\t mkdir ${workdir}/${MLmethod}/asymPlots"
mkdir -p $workdir/$MLmethod/asymPlots
printgreen "\t mkdir ${workdir}/${MLmethod}/binnedPlots"
mkdir -p $workdir/$MLmethod/binnedPlots
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
    printgreen "\t Generating /volatile link in $workdir/$MLmethod"
    ln -s $volatiledir/$MLmethod $workdir/$MLmethod/linkToVolatile
    printgreen "\t Generating MLinputs.yaml file in $workdir/$MLmethod"
    cp $NNPI0DIR/catboost/input/input_noresonance.yaml $workdir/$MLmethod/MLinputs.yaml
    # Create the train.sh script
    printgreen "\t Generating train.sh file in $workdir/$MLmethod"
    cat >> $workdir/$MLmethod/train.sh << EOF
cd $NNPI0DIR/catboost/
python train.py --data_dir ${volatiledir}/${MLmethod}/MLinput \
--depth 3 \
--make_plots true \
--subdata MC_inbending \
--model_dir $CATBOOST_MODEL_INBENDING \
--input_yaml $CATBOOST_INPUT_YAML

python train.py --data_dir ${volatiledir}/${MLmethod}/MLinput \
--depth 3 \
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
    --preprocess_data_dir ${volatiledir}/${MLmethod}/MLinput \
    --subdata RGA_inbending \
    --output_dir ${volatiledir}/${MLmethod}/MLoutput
    
python predict.py \
    --version nSidis \
    --model_dir $CATBOOST_MODEL_OUTBENDING \
    --preprocess_data_dir ${volatiledir}/${MLmethod}/MLinput \
    --subdata RGA_outbending \
    --output_dir ${volatiledir}/${MLmethod}/MLoutput
    
python predict.py \
    --version MC \
    --model_dir $CATBOOST_MODEL_INBENDING \
    --preprocess_data_dir ${volatiledir}/${MLmethod}/MLinput \
    --subdata MC_inbending \
    --output_dir ${volatiledir}/${MLmethod}/MLoutput
    
python predict.py \
    --version MC \
    --model_dir $CATBOOST_MODEL_OUTBENDING \
    --preprocess_data_dir ${volatiledir}/${MLmethod}/MLinput \
    --subdata MC_outbending \
    --output_dir ${volatiledir}/${MLmethod}/MLoutput
    
cd $pwd
EOF

    # Create runSidis.sh script
    printgreen "\t Generating runSidis.sh file in $workdir/$MLmethod"
    cat >> $workdir/$MLmethod/runSidis.sh << EOF
#!/bin/bash
files=($(basename $SCIPIODIR/macros/process/* ))
script=\$1
if [ \$# -lt 1 ]; then
    echo "./runSidis.sh <PROCESS SCRIPT>"
    echo -e "\n Please use one of the available files\n \t\$files"
    exit 1
fi

if [ ! -f "$SCIPIODIR/macros/process/\$script" ]; then
    echo "Error: script not found."
    echo -e "\n Please use one of the available files\n \t\$files"
    exit 1
fi

cd $SCIPIODIR/macros/slurm

chmod +x runSidis.sh

./runSidis.sh nSidis catboost ${volatiledir}/${MLmethod}/MLoutput ${volatiledir}/${MLmethod}/postprocess \$script
./runSidis.sh MC catboost ${volatiledir}/${MLmethod}/MLoutput ${volatiledir}/${MLmethod}/postprocess \$script

cd $pwd
EOF
    
    # Create runBin.sh script
    printgreen "\t Generating Binning.yaml file in $workdir/$MLmethod"
    cp $SCIPIODIR/utils/Binning.yaml $workdir/$MLmethod/Binning.yaml
    printgreen "\t Generating runBin.sh file in $workdir/$MLmethod"
    cat >> $workdir/$MLmethod/runBin.sh << EOFmain
#!/bin/bash
num=\$(grep -o "parentDirectory" Binning.yaml | wc -l)
v=("MC" "nSidis")
vidx=(1 0)
for ((i=0; i < \$num; i++))
do
    for ((j=0; j < \${#vidx[@]}; j++))
    do
        version=\${v[\$j]}
        isMC=\${vidx[\$j]}
        
        slurmshell=${slurmdir}"/binning_\${i}_\${j}.sh"
        
        slurmslurm=${slurmdir}"/binning_\${i}_\${j}.slurm"

        touch \$slurmshell
        touch \$slurmslurm
        chmod +x \$slurmshell

        cat >> \$slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=job_binning_$version_\${i}_\${j}
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=${logdir}/binning_\${i}_\${j}.out
#SBATCH --error=${logdir}/binning_\${i}_\${j}.err
$slurmshell
EOF
        cat >> \$slurmshell << EOF
source /group/clas12/packages/setup.csh
module load clas12/pro
clas12root -q '/work/clas12/users/gmat/scipio/macros/analysis/binner.C("${volatiledir}/${MLmethod}/postprocess/\$version*.root","${volatiledir}/${MLmethod}/postprocess_binned","${workdir}/${MLmethod}/Binning.yaml",\$isMC,\$i)'
EOF
    
    sbatch \$slurmslurm
    done
done
EOFmain

    # Create runBru.sh script
    printgreen "\t Generating runBru.sh file in $workdir/$MLmethod"
    cat >> $workdir/$MLmethod/runBru.sh << EOFmain
#!/bin/bash                                                                                                                                                                  
input_dir=${volatiledir}/${MLmethod}/postprocess_binned
L=2
threshold=0.9
Mggmin=0.07
Mggmax=0.22
sidebandMin=0.22
sidebandMax=0.4
BRUFIT=$SCIPIODIR/brufit


for dir in \$(find "\$input_dir" -mindepth 1 -maxdepth 1 -type d)
do
  for file in \$(find "\$dir" -name "*.root")
  do
    dirname=\$(basename \$dir)
    filename=\$(basename \$file)
    isMC=0
    if [[ \$filename == MC* ]]; then
        isMC=1
    fi
    if [[ ! \$filename =~ ^(MC|nSidis).* ]]; then
        continue
    fi
    programOutput=\$(clas12root -b -q -l /work/clas12/users/gmat/scipio/src/ReadTTrees.C\\(\\"\$file\\"\\))
    ttrees=\${programOutput#*...}
    ttrees=\${ttrees#?}
    
    echo "Submitting job for TTrees in \$filename"
    echo -e "\t \$ttrees"
    echo -e "\n"
    
    # Create slurm files
    for ttree in \${ttrees[@]}
    do
        if [[ \$isMC == 0 ]]; then
            versions=("all" "Fall2018_inbending" "Fall2018_outbending" "Spring2019_inbending")
        else
            versions=("all" "MC_inbending" "MC_outbending")
        fi
        for version in \${versions[@]}
        do
            slurmshell=${slurmdir}"/\${dirname}_\${filename}_\${isMC}_\${ttree}_\${version}.sh"
            slurmslurm=${slurmdir}"/\${dirname}_\${filename}_\${isMC}_\${ttree}_\${version}.slurm"

            touch \$slurmshell
            touch \$slurmslurm
            chmod +x \$slurmshell

            cat >> \$slurmslurm << EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=${memPerCPU}
#SBATCH --job-name=job_\${dirname}_\${filename}_\${isMC}_\${ttree}_\${version}
#SBATCH --cpus-per-task=${nCPUs}
#SBATCH --time=24:00:00
#SBATCH --output=${logdir}/\${dirname}_\${filename}_\${isMC}_\${ttree}_\${version}.out
#SBATCH --error=${logdir}/\${dirname}_\${filename}_\${isMC}_\${ttree}_\${version}.err
\$slurmshell
EOF

            cat >> \$slurmshell << EOF
#!/bin/tcsh
source /group/clas12/packages/setup.csh
module load clas12/pro
/u/site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/bin/root \$BRUFIT/macros/LoadBru.C -b -q -l $SCIPIODIR/macros/analysis/bruana_pipluspi0.C\\(\\"\$input_dir/\$dirname\\",\\"\$file\\",\\"\$ttree\\",\$L,\$threshold,\$Mggmin,\$Mggmax,\$sidebandMin,\$sidebandMax,\$isMC,\\"\$version\\"\\)
/u/site/12gev_phys/2.4/Linux_CentOS7.7.1908-gcc9.2.0/root/6.20.04/bin/root \$BRUFIT/macros/LoadBru.C -b -q -l $SCIPIODIR/macros/analysis/bruana_pipluspi0.C\\(\\"\$input_dir/\$dirname\\",\\"\$file\\",\\"\$ttree\\",\$L,1,\$Mggmin,\$Mggmax,\$sidebandMin,\$sidebandMax,\$isMC,\\"\$version\\"\\)
EOF

            sbatch \$slurmslurm
        done
    done
done
done
EOFmain
fi

#######################################################################################################
#######################################################################################################
#######################################################################################################


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
