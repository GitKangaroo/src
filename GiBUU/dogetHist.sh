export GFSa3nuAr="outplot/outAna3_GiBUUT2KGFSa3nuAr_FinalEventsList_005_T2K_nu_freeDelta_T0_noFSI_Flux9_noEuCut_Argon_1kruns.root"
export GFSa0nuC="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_MINERvAGiBUUGFSa0nuC_FinalEventsList_005_MINERvA_nu_freeDelta_T0_FSI_Flux25_noEuCut_1kruns_14_10jobs.root"
export GFSa0nuC_lu="/data/t2k/xlu/.analysis/anaGiBUU/outplot/outAna0_nuCarbonMINERvAT0GiBUU2019a0newbkgparticle_FinalEventsList_005_MINERvA_nu_freeDelta_T0_FSI_Flux25_noEuCut_1kruns_2_998jobs.root"
export GFSa0nuAr="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUGFSa0nuAr_rep_FinalEventsList_005_MINERvA_nu_freeDelta_T0_FSI_Flux25_noEuCut_Argon_1kruns_13_10jobs.root"
export GFSa0nuC="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUGFSa0nuC_rep_FinalEventsList_005_MINERvA_nu_freeDelta_T0_FSI_Flux25_noEuCut_1kruns_14_10jobs.root"

export GFSa0nuC2019="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_MINERvAGiBUUGFSa0nuC2019_FinalEvents_1000runs_911jobs_MINERvA-nu_freeDelta_T0_FSI_Flux25_noEuCut_1kruns.root"

export GFSa0nuC_DUNE="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUDUNEGFSa0nuC_DUNE_FinalEventsList_005_DUNE_nu_freeDelta_T0_FSI_Flux15_noEuCut_Carbon_1kruns_1_999jobs.root"
export GFSa0nuAr_DUNE="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUDUNEGFSa0nuAr_DUNE_FinalEventsList_005_DUNE_nu_freeDelta_T2_FSI_Flux15_noEuCut_Argon_1kruns_2_997jobs.root"

export Step1_Carbon_MINERvA_ens40_nrun1000_T0="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUStep1_Carbon_MINERvA_ens40_nrun1000_T0_Step1.root"
export Step2_Carbon_MINERvA_ens4000_nrun1_T0="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUStep2_Carbon_MINERvA_ens4000_nrun1_T0_Step2.root"
export Step3_Carbon_DUNE_ens4000_nrun1_T0="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUStep3_Carbon_DUNE_ens4000_nrun1_T0_Step3.root"
export Step4_Argon_MINERvA_ens40_nrun1000_T0="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUStep4_Argon_MINERvA_ens40_nrun1000_T0_Step4.root"
export Step5_Argon_MINERvA_ens4000_nrun1_T0="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUStep5_Argon_MINERvA_ens4000_nrun1_T0_Step5.root"
export Step6_Argon_DUNE_ens4000_nrun1_T0="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUStep6_Argon_DUNE_ens4000_nrun1_T0_Step6.root"
export Step7_Argon_DUNE_ens4000_nrun1_T2="/data/t2k/yangk/analysis/anaGiBUU/outplot/outAna0_GiBUUStep7_Argon_DUNE_ens4000_nrun1_T2_Step7.root"

tag=GiBUU
./exe_mk.sh getHist
echo exe_mk.sh done

#opt=${tag}GFSa0nuC_Hist_CC0piNpID;
#./getHist $GFSa0nuC          ${opt} 1 > see${opt}.log &

#opt=${tag}GFSa0nuAr_Hist_CC0piNpID;
#./getHist $GFSa0nuAr          ${opt} 1 > see${opt}.log &

#GFS 1 CC0piNpID
opt=${tag}GFSa0nuC_Hist_CC0piNpID_DUNE;
./getHist $GFSa0nuC_DUNE          ${opt} 1 > see${opt}.log &

#GFS 1 CC0piNpID
opt=${tag}GFSa0nuAr_Hist_CC0piNpID_DUNE;
./getHist $GFSa0nuAr_DUNE          ${opt} 1 > see${opt}.log &


exit
opt=${tag}Step1_Carbon_MINERvA_ens40_nrun1000_T0_Hist;
./getHist $Step1_Carbon_MINERvA_ens40_nrun1000_T0          ${opt} 1 > see${opt}.log &

opt=${tag}Step2_Carbon_MINERvA_ens4000_nrun1_T0_Hist;
./getHist $Step2_Carbon_MINERvA_ens4000_nrun1_T0          ${opt} 1 > see${opt}.log &

opt=${tag}Step3_Carbon_DUNE_ens4000_nrun1_T0_Hist;
./getHist $Step3_Carbon_DUNE_ens4000_nrun1_T0          ${opt} 1 > see${opt}.log &

opt=${tag}Step4_Argon_MINERvA_ens40_nrun1000_T0_Hist;
./getHist $Step4_Argon_MINERvA_ens40_nrun1000_T0          ${opt} 1 > see${opt}.log &

opt=${tag}Step5_Argon_MINERvA_ens4000_nrun1_T0_Hist;
./getHist $Step5_Argon_MINERvA_ens4000_nrun1_T0          ${opt} 1 > see${opt}.log &

opt=${tag}Step6_Argon_DUNE_ens4000_nrun1_T0_Hist;
./getHist $Step6_Argon_DUNE_ens4000_nrun1_T0          ${opt} 1 > see${opt}.log &

opt=${tag}Step7_Argon_DUNE_ens4000_nrun1_T2_Hist;
./getHist $Step7_Argon_DUNE_ens4000_nrun1_T2          ${opt} 1 > see${opt}.log &

exit

#GFS 1 CC0piNpID
opt=${tag}GFSa0nuC_Hist_CC0piNpID_DUNE;
./getHist $GFSa0nuC_DUNE          ${opt} 1 > see${opt}.log &
exit
#GFS 3 CC1piNpID
opt=${tag}GFSa0nuC_Hist_CC1piNpID;
./getHist $GFSa0nuC          ${opt} 2 > see${opt}.log &

#GFS 3 LOWRECOIL
opt=${tag}GFSa0nuC_Hist_LOWRECOIL;
./getHist $GFSa0nuC          ${opt} 8 > see${opt}.log &

#GFS 3 LOWRECOIL0piNp
opt=${tag}GFSa0nuC_Hist_LOWRECOIL0piNp;
./getHist $GFSa0nuC          ${opt} 9 > see${opt}.log &

#GFS 3 NUBAR1PI
opt=${tag}GFSa0nuC_Hist_NUBAR1PI;
./getHist $GFSa0nuC          ${opt} 16 > see${opt}.log &

#GFS 3 MMECCQE
opt=${tag}GFSa0nuC_Hist_MMECCQE;
./getHist $GFSa0nuC          ${opt} 32 > see${opt}.log &












