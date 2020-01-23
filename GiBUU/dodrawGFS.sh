#!/bin/bash

exe=drawGFS; 


mkdir outplot/final -p

root -b -l <<EOF
.L style.cxx+
.L ${exe}.C++
.q
EOF

     
for tag in "GiBUUC"
#"GENIEOOB T2K" "GENIEOOB MINERvA"
#"GENIE MINERvA"
#"GiBUU MME"
#
#"GiBUUnoFSI T2K"
#"GiBUU T2K" "NuWro T2K" "GiBUU MINERvA" #final_GFS1piT2KGiBUUnu_additional_GiBUU_MINERvA_vs_T2K
#"GiBUU MINERvA"
#"GiBUU MINERvA" "NuWro T2K" 
#no use, already set as other "NuWro MINERvA"
#"GiBUU MINERvA" 
#
#"NuWro T2K" "GiBUU MINERvA"
#
do
echo tag $tag

for mode in 0
#0 2
#2 
#5 6 11
#2
#$(seq 11 14)
#5 6 
#$(seq 0 4)
#5 6 15 16
#
#4
#
#15
#13 14
#9 10
#$(seq 5 8)
#7 8
#5 6
#12
#
#$(seq 5 8)
#$(seq 7 8)
#8
#$(seq 0 7)
do
    
    for ii in $(seq 0 37 )
#$(seq 0 34)
    do
        
        if [ ${mode} -lt 2 ]
        then
            if [ ${ii} -gt 7 ]
            then
                continue
            fi
        elif [ ${mode} -lt 5 ]
        then
            if [ ${ii} -gt 13  -a ${ii} -lt 35  ]
            then
                continue
            fi
        elif [ ${mode} -lt 9 -o ${mode} = 15 -o ${mode} = 16 ]
        then
            if [ ${ii} -lt 14 -o ${ii} -gt 20 ]
            then
                continue
            fi
        elif [ ${mode} -lt 11 ]
        then
            if [ ${ii} -lt 15 -o ${ii} -gt 26 ]
            then
                continue
            fi
        elif [ $mode -lt 13 ]
        then
            if [ ${ii} -lt 28 ]
            then 
                continue
            fi
        elif [ $mode -lt 15 ]
        then
            if [ ${ii} != 27 ]
            then 
                continue
            fi
        fi
                

        echo mode $mode ii $ii
        
        root -b -l <<EOF
.L style.cxx+
.L ${exe}.C+
${exe}(${mode}, $ii, "${tag}")
.q
EOF
        
    done
done

done
