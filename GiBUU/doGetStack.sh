#remember to change var correspondingly

echo
echo "check normalization: GFS1pi is carbon only normalization (1/12), others should be 1/13"
echo
anaid=GFS0piDUNEGiBUUC
#anaid=GFS0piMINERvAGiBUUC_re
#MMECCQEMINERvAMEGiBUU
#GFSPIZEROMINERvAGiBUU
#GFSEXP1piT2KGiBUU
#
#GFS0piMMEGiBUU
#NUBAR1PIMINERvAGiBUU
#LOWRECOILMINERvAGiBUU
#GFS1piMINERvAGiBUU
#GFS0piMINERvAGiBUU
#
#
#GFSEXP1piT2KGiBUUnoFSI
#GFSEXP1piT2KGiBUU
#
#GFS1piT2KGiBUU
#
#GFS0piT2KGiBUU

#GFS1piT2KGiBUUnoFSI
#GFS1piT2KGiBUU
#NUBAR1PIMINERvAGiBUU
#NUBAR1PIsigmaEnuMINERvAGiBUU
#
#GFS0piMINERvAGENIE
#
#
#GFS0piT2KGENIE
#
#GFS0piT2KGiBUU
#GFS0piMINERvAGiBUU
#
#$1
#LOWRECOILCutR_GiBUUMINERvA
#GFS1piMINERvAGiBUU
#LOWRECOILCugnonFSI_GiBUUMINERvA
#GFS1piT2KNuWro
#GFS0piT2KNuWro
#
#GFS1piMINERvANuWro
#GFS0piMINERvANuWro
#GFS0piMINERvAGiBUU
#
#GFS1pi
#
#NUBAR1PIsigmaEnu
#
#GFS0piT2KT1
#GFS0piT2K
#GFS0piT1
#
#LOWRECOIL
#LOWRECOILnoFSI
#LOWRECOIL0piNp
#
#NUBAR1PI
#NUBAR1PInoFSI
#NUBAR1PIsigmaEnu
#NUBAR1PIsigmaEnunoFSI

mkdir outStack/${anaid} -p
ls outStack/${anaid}
echo
echo

for locdir in ${anaid}nu{,bar} 
#${anaid}nu_T0                                                       
#{nubar,nubar_noFSI}                                                 
#${anaid}nu                                                          
#${anaid}nu_T{0,1}                                                   
#${anaid}_NuWroMINERvAnu{,bar} ${anaid}_GiBUUMINERvAnu{,bar}         
#${anaid}_GiBUUMINERvAnu                                             
#${anaid}_GiBUUMINERvAnubar                                          
#${anaid}_GiBUUMINERvAnu{,bar}  
do
    echo locdir $locdir
    echo
    fin=outplot/${locdir}.root

    if [ ! -e $fin ]
    then
        continue
    fi

    ls $fin
    echo
    tag=$(basename $fin | awk -F\. '{print $1}')
    echo tag $tag
    echo

    commonVar="muonmomentum muontheta enu Q2 xBj xrest Wtrue Wrest "
    if [ $anaid == GFS0piDUNEGiBUUC -o $anaid == GFS0piMINERvAGiBUU -o $anaid == GFS0piMINERvAGiBUUC_re -o $anaid == GFS1piMINERvAGiBUU -o $anaid == GFS0piMMEGiBUU -o $anaid == GFSEXP1piT2KGiBUU -o $anaid == GFSPIZEROMINERvAGiBUU ]
    then
# dpTT dpTy
        varray=${commonVar}"dphit protonmomentum protontheta dalphat dpt neutronmomentum pionEk pionmomentum piontheta baryonmomentum baryontheta baryonmass dpTT protonTT pionTT "
    elif [ $anaid == LOWRECOILMINERvAGiBUU ]
    then
        varray=${commonVar}$(echo eav{0,1,2,3,4,5} q3)
    elif [ $anaid == MINERvAMEGiBUUMMECCQE ]
#MMECCQEMINERvAMEGiBUU ]
    then
        varray=${commonVar}$(echo  q2qe muonpt mupz muonptVSmupz mupzVSmuonpt)
#muonpt{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14} mupz{0,1,2,3,4,5,6,7,8,9,10,11,12}
    elif [ $anaid == NUBAR1PIMINERvAGiBUU ]
    then
        varray=${commonVar}"pionEk piontheta pionmomentum"
    else
        echo wrong anaid $anaid
        exit
    fi

    echo $varray

    #continue

    for var in $varray
    do

        if [ $var = neutronmomentum ]; then nuwrovar=pn; fi
        if [ $var = dphit ]; then nuwrovar=phiT; fi
        if [ $var = muonmomentum ]; then nuwrovar=muonmom; fi
        if [ $var = muontheta ]; then nuwrovar=muonangle; fi
        if [ $var = protonmomentum ]; then nuwrovar=nucmom; fi
        if [ $var = protontheta ]; then nuwrovar=nucangle; fi
        if [ $var = dalphat ]; then nuwrovar=alphaT; fi
        if [ $var = dpt ]; then nuwrovar=dpT; fi

        if [ $var = pionEk ]; then nuwrovar=pionenkin; fi
        if [ $var = pionmomentum ]; then nuwrovar=pionmom; fi
        if [ $var = piontheta ]; then nuwrovar=pionangle; fi
        if [ $var = baryonmomentum ]; then nuwrovar=piNmom; fi
        if [ $var = baryontheta ]; then nuwrovar=piNangle; fi
        if [ $var = baryonmass ]; then nuwrovar=invmass; fi

        echo var $var nuwrovar $nuwrovar
        echo
    
        root -b -l -q 'getStack.C("'${fin}'", "'${anaid}'","'${tag}'","'${var}'","'${nuwrovar}'")'

    done
done

cd outStack/${anaid}
ls *.root
hadd ${anaid}.root *.root > seemerge.log


