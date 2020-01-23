#anaid_Ar=GFS0piMINERvAGiBUUAr_re
#anaid_C=GFS0piMINERvAGiBUUC_re
anaid_Ar=GFS0piDUNEGiBUUAr
anaid_C=GFS0piDUNEGiBUUC

tag=GFS0piDUNEGiBUU

mkdir outStackSame
ls outStackSame
echo
echo

fin_Ar=outStack/${anaid_Ar}/${anaid_Ar}.root
fin_C=outStack/${anaid_C}/${anaid_C}.root

commonVar="muonmomentum muontheta enu Q2 xBj xrest Wtrue Wrest "

varray=${commonVar}"dphit protonmomentum protontheta dalphat dpt neutronmomentum pionEk pionmomentum piontheta baryonmomentum baryontheta baryonmass dpTT protonTT pionTT "

for var in $varray
do

   echo var $var
   echo

   root -b -l -q 'drawSame.C("'${fin_Ar}'", "'${fin_C}'","'${tag}'","'${var}'")'
done

