dir0=$(pwd)
jobid=1
nrepeat=5
#125

#aftertag="CugnonFSI/"

for rawjobcard in noEuCut/005_MINERvA_nu_freeDelta_T0_FSI_Flux25_noEuCut_Argon_1kruns.job
#for rawjobcard in noEuCut/005_T2K_nu_freeDelta_T0_noFSI_Flux9_noEuCut_hydrogen_1kruns.job
#noEuCut/{005_T2K_nu_freeDelta_T0_noFSI_Flux9_noEuCut_hydrogen_1kruns.job,005_T2K_nu_freeDelta_T0_noFSInoFSInoFSInoFSI_Flux9_noEuCut_1kruns.job}
#sigmaEnu/*/005*.job
#EuCut1500to10000/005_MINERvA_nubar_freeDelta_T0_*
#EuCut2to6/005_MINERvA_nu_freeDelta_T0_noFSI_Flux25_EuCut2to6_hydrogen_1kruns.job
#20190613 -> EuCut2to6/005_MINERvA_nu*.job
#ECT2019June/012_HiPion_N.job
#012_HiPion_A.job
#noEuCut/005*T2K*.job
do
    
batchtag=$(basename $rawjobcard | awk -F\. '{print $1}')
echo $batchtag

dirtag=$(dirname $rawjobcard )/output/Step4/${aftertag}test_${batchtag}_${jobid}
echo $dirtag

jobcard=$(readlink -f ${rawjobcard})
echo ${jobcard}

gibuuexe=$(readlink -f GiBUU.x)
echo GiBUU exe $gibuuexe
echo GiBUU jobcard $jobcard

for kk in $(seq 1 ${nrepeat})
do
    dir=${dir0}/${dirtag}_out${kk}

    echo ${dir}
    
    mkdir -p ${dir}
        
    cd ${dir}

    cp ${jobcard} .

    jobfile=s${jobid}.sh
###
    cat > ${jobfile} <<EOF

cd ${dir}
    source /home/yangk/.profile
    job1="${gibuuexe} < $(basename ${jobcard}) > see_${batchtag}.log"
    echo \$job1
date

    eval \$job1

date
EOF
##
    chmod +x ${jobfile}
    #qsub -lnodes=1:ppn=3 -lcput=00:10:00 -N $jobfile  ./$jobfile 
    #qsub -lcput=24:00:00  -N $jobfile  ./$jobfile  #10k runs take 12.5h
    #qsub -lcput=00:40:00  -N $jobfile  ./$jobfile  #10k runs take 12.5h
#120TimeStep
    #qsub -lcput=3:20:00  -N $jobfile  ./$jobfile  #10k runs take 12.5h
#150TimeStep
    qsub -lcput=3:10:00  -N $jobfile  ./$jobfile  #10k runs take 12.5h
    #qsub -lcput=72:00:00  -N $jobfile  ./$jobfile 
    #nohup ./$jobfile &
    #./$jobfile

    cd $dir0
    
    jobid=$( expr ${jobid} + 1 )

    echo
done

done

