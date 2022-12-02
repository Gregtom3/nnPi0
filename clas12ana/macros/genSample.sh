#!/bin/bash
hl="--------------------------------------------------------------------------------------------------------------------------------------------"
nFiles=2
nEvents=10000 # per file
ana="MC" # either MC or SIDIS
preprocess="catboost" # catboost or particleNet
volatiledir=/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML



if [ $ana == "MC" ]; then
    declare -a hipofiles=("/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg*/*.hipo" "/cache/clas12/rg-a/production/montecarlo/clasdis/fall\2018/torus+1/v1/bkg*/*.hipo")
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

    echo $hl
    echo "Reading hipo file $hipo"
    echo $hl
    
    if [ -f "$raw_out" ]; then
	echo "$raw_out exists...skipping..."
    else
	clas12root -b -q "pi0_readHipo.C(\"${hipo}\",\"${raw_out}\",$beamE,$nEvents,1)"
	echo $hl
    fi

    echo "Preprocessing hipo file $hipo"
    echo $hl
    clas12root -b -q "pi0_preprocess_${preprocess}.C(\"${raw_out}\",\"${preprocess_out}\")"
    echo $hl
    
    filenum=$((filenum+1))
    if [ $filenum -gt $nFiles ]; then
	echo $hl
	break
    fi
done
done
