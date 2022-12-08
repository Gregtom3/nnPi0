#!/bin/bash
USERNAME="$USER"
hl="--------------------------------------------------------------------------------------------------------------------------------------------"
nFiles=5
nEvents=1000000 # per file
ana="nSidis" # either MC or nSidis
preprocess="catboost" # catboost or particleNet
volatiledir=/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML

slurm=true
logdir=""
slurmdir=""
if [ $slurm ]; then 
    
        NOW=$( date '+%F_%H_%M' )
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
    preprocess_out=$volatiledir"/preprocess_${preprocess}/${ana}_${runNumber}.root"
    

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
        
        cat >> $slurmshell << EOF
#!/bin/tcsh
source /group/clas12/packages/setup.csh
module load clas12/pro
clas12root -b -q "pi0_readHipo.C(\"${hipo}\",\"${raw_out}\",$beamE,$nEvents,$hipo_is_mc)"
clas12root -b -q "pi0_preprocess_${preprocess}.C(\"${raw_out}\",\"${preprocess_out}\")"
echo "Done"
EOF
        
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
