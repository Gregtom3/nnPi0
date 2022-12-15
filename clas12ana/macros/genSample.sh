#!/bin/bash
USERNAME="$USER"
hl="--------------------------------------------------------------------------------------------------------------------------------------------"

if [ $# -lt 2 ]; then
  echo """
  USAGE: $0 [MC/nSidis] [ML Method] [flags(optional)]
  Automates the generation of raw TTrees from the CLAS12 HIPO data
  After the raw TTrees are created (or if they are created already) this program preprocesses
  the raw TTrees into the input space of various ML algorithms
   - [MC/nSidis]: Specifies if we are analyzing MC or nSidis data
   - [ML Method]: Options (catboost particleNet)            
   - [flags]:
              -o      (overwrite existing raw TTrees)
              -nFiles <INT> (maximum number of files per hipo directory: default=5)
              -noslurm (run all the processing code in command line without slurm)
              -maxEvents <INT> (set the maximum number of events analyzed per hipo file)
  """
  exit 2
fi

PREPROCESSES="catboost particleNet"
ANAS="MC nSidis"
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

ana=${args[0]}
preprocess=${args[1]}
nFiles=5

if [ ! -z ${flags["nFiles"]} ]; then
    nFiles=${flags["nFiles"]}
fi

if ! echo "$PREPROCESSES" | grep -q "$preprocess"; then
    echo $hl
    echo "ERROR: ML Method $preprocess not valid. Must use from following list [$PREPROCESSES]"
    exit 2
else
    rungroup="${rungroup//-}"
fi

if ! echo "$ANAS" | grep -q "$ana"; then
    echo $hl
    echo "ERROR: Analysis type $ana not valid. Must use from following list [$ANAS]"
    exit 2
fi

nEvents=1000000 # per file
if [ ! -z ${flags["maxEvents"]} ]; then
    nEvents=${flags["maxEvents"]}
fi

volatiledir=/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML
pwd=$PWD

slurm=true
if [ ! -z ${booleans["noslurm"]} ]; then
    slurm=false
fi



logdir=""
slurmdir=""
if [ $slurm ]; then 
    
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

    raw_out=$volatiledir"/raw/${ana}_${runNumber}.root"
    preprocess_out=$volatiledir"/$preprocess/preprocess_pi0/${ana}_${runNumber}.root"

    if [ $slurm ]; then 
    
        slurmshell=${slurmdir}"/${ana}_${runNumber}_${preprocess}.sh"
        slurmslurm=${slurmdir}"/${ana}_${runNumber}_${preprocess}.slurm"
            
        touch $slurmshell
        touch $slurmslurm
        chmod +x $slurmshell
        
        cat >> $slurmslurm <<EOF
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=${memPerCPU}
#SBATCH --job-name=job_${ana}_${runNumber}_${preprocess}
#SBATCH --cpus-per-task=${nCPUs}
#SBATCH --time=24:00:00
#SBATCH --output=${logdir}/${ana}_${runNumber}_${preprocess}.out
#SBATCH --error=${logdir}/${ana}_${runNumber}_${preprocess}.err
$slurmshell
EOF

        if [ ! -z ${booleans["o"]} ]; then
            cat >> $slurmshell << EOF
#!/bin/tcsh
source /group/clas12/packages/setup.csh
module load clas12/pro
clas12root -b -q $PWD/pi0_readHipo.C\(\"${hipo}\",\"${raw_out}\",$beamE,$nEvents,$hipo_is_mc\)
clas12root -b -q $PWD/pi0_preprocess_${preprocess}.C\(\"${raw_out}\",\"${preprocess_out}\"\)
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
clas12root -b -q $PWD/pi0_preprocess_${preprocess}.C\(\"${raw_out}\",\"${preprocess_out}\"\)
echo "Done"
EOF
	fi

        sbatch $slurmslurm
        
    else    

        echo $hl
        echo "Reading hipo file $hipo"
        echo $hl

        if [ -f "$raw_out" ]; then
        echo "$raw_out exists...skipping..."
        else
        clas12root -b -q "pi0_readHipo.C(\"${hipo}\",\"${raw_out}\",$beamE,$nEvents,$hipo_is_mc)"
        fi

        echo "Preprocessing hipo file $hipo"
        echo $hl
        clas12root -b -q "pi0_preprocess_${preprocess}.C(\"${raw_out}\",\"${preprocess_out}\")"
        echo $hl

    fi
    
    filenum=$((filenum+1))
    if [ $filenum -eq $nFiles ]; then
        break
    fi
    
done
done


echo $hl
echo "Location of log files: $log_dir"
echo "Location of slurm files: $slurm_dir"
