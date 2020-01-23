for ii in noEuCut/output/tes*
do

    jj=$(head -n 1 $ii/FinalEvents.dat)bb

    if [ "$jj" == "bb" ]
    then
        echo empty $ii
        mv $ii bad/
    fi
done

for ii in noEuCut/output/tes*
do 
    #jj=$(grep "BUU simulation: finished" $ii/see*.log)bb 

    jj=$(tail $ii/see*.log -n 5 | head -n 1 | awk '{if($4=="finished")print}')bb

    if [ "$jj" == "bb" ]
    then
        echo bad $ii
        mv $ii bad/
    fi

done

echo done
