#!/bin/sh

unalias ls

path_of_this_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
HOMEBASE=${path_of_this_dir}

fake=0
mt=1

#mupath=/group/atlas/data/tauD3PD/zprimetautau/lephad/2014_muhad_paper_hists/2014-08-05
#epath=/group/atlas/data/tauD3PD/zprimetautau/lephad/2014_ehad_paper_hists/2014-08-13
mupath=/group/atlas/data/tauD3PD/zprimetautau/lephad/2014_muhad_paper_hists/2014-11-03
epath=/group/atlas/data/tauD3PD/zprimetautau/lephad/2014_ehad_paper_hists/2014-10-22

if [ $mt == 1 ] ; then
    mumainpath=`ls -d ${mupath}/*MAIN`
    mumuonscaleuppath=`ls -d ${mupath}/*MUON_SCALE_UP`
    mumuonscalednpath=`ls -d ${mupath}/*MUON_SCALE_DN`
    mumuoniduppath=`ls -d ${mupath}/*MUON_ID_UP`
    mumuonmsuppath=`ls -d ${mupath}/*MUON_MS_UP`
    mutesrealuppath=`ls -d ${mupath}/*TESREAL_UP`
    mutesreakdnpath=`ls -d ${mupath}/*TESREAL_DN`
    mutesfakeuppath=`ls -d ${mupath}/*TESFAKE_UP`
    mutesfakednpath=`ls -d ${mupath}/*TESFAKE_DN`
    mujesuppath=`ls -d ${mupath}/*JES_UP`
    mujesdnpath=`ls -d ${mupath}/*JES_DN`
    mumetscalesofttermsuppath=`ls -d ${mupath}/*MET_ScaleSoftTerms_UP`
    mumetscalesofttermsdnpath=`ls -d ${mupath}/*MET_ScaleSoftTerms_DN`
    mumetresosofttermsuppath=`ls -d ${mupath}/*MET_ResoSoftTerms_UP`
    mumueffuppath=`ls -d ${mupath}/*MU_EFF_UP`
    mumueffdnpath=`ls -d ${mupath}/*MU_EFF_DN`
    mutau3puppath=`ls -d ${mupath}/*TAU3P_UP`
    mutau3pdnpath=`ls -d ${mupath}/*TAU3P_DN`
    mumutriguppath=`ls -d ${mupath}/*MU_TRIG_UP`
    mumutrigdnpath=`ls -d ${mupath}/*MU_TRIG_DN`
    mutauismuonuppath=`ls -d ${mupath}/*TAU_ISMUON_UP`
    mutauismuondnpath=`ls -d ${mupath}/*TAU_ISMUON_DN`
    mutauiduppath=`ls -d ${mupath}/*TAUID_UP`
    mutauiddnpath=`ls -d ${mupath}/*TAUID_DN`
    mukfactoruppath=`ls -d ${mupath}/*KFACTOR_UP`
    mukfactordnpath=`ls -d ${mupath}/*KFACTOR_DN`
    emainpath=`ls -d ${epath}/*MAIN`
    eelescaleuppath=`ls -d ${epath}/*ELE_SCALE_UP`
    eelescalednpath=`ls -d ${epath}/*ELE_SCALE_DN`
    eeleresuppath=`ls -d ${epath}/*ELE_RES_UP`
    etesrealuppath=`ls -d ${epath}/*TESREAL_UP`
    etesreakdnpath=`ls -d ${epath}/*TESREAL_DN`
    etesfakeuppath=`ls -d ${epath}/*TESFAKE_UP`
    etesfakednpath=`ls -d ${epath}/*TESFAKE_DN`
    ejesuppath=`ls -d ${epath}/*JES_UP`
    ejesdnpath=`ls -d ${epath}/*JES_DN`
    emetscalesofttermsuppath=`ls -d ${epath}/*MET_ScaleSoftTerms_UP`
    emetscalesofttermsdnpath=`ls -d ${epath}/*MET_ScaleSoftTerms_DN`
    emetresosofttermsuppath=`ls -d ${epath}/*MET_ResoSoftTerms_UP`
    eeleeffuppath=`ls -d ${epath}/*ELE_EFF_UP`
    eeleeffdnpath=`ls -d ${epath}/*ELE_EFF_DN`
    etau3puppath=`ls -d ${epath}/*TAU3P_UP`
    etau3pdnpath=`ls -d ${epath}/*TAU3P_DN`
    eeletriguppath=`ls -d ${epath}/*ELE_TRIG_UP`
    eeletrigdnpath=`ls -d ${epath}/*ELE_TRIG_DN`
    etauiselectronuppath=`ls -d ${epath}/*TAU_ISELECTRON_UP`
    etauiselectrondnpath=`ls -d ${epath}/*TAU_ISELECTRON_DN`
    etauiduppath=`ls -d ${epath}/*TAUID_UP`
    etauiddnpath=`ls -d ${epath}/*TAUID_DN`
    ekfactoruppath=`ls -d ${epath}/*KFACTOR_UP`
    ekfactordnpath=`ls -d ${epath}/*KFACTOR_DN`

    paper_input_mttot.py --mufolder=${mupath} --mu_main=${mumainpath} --mu_muonscaleup=${mumuonscaleuppath} --mu_muonscaledn=${mumuonscalednpath} --mu_muonidup=${mumuoniduppath} --mu_muonmsup=${mumuonmsuppath} --mu_tesrealup=${mutesrealuppath} --mu_tesrealdn=${mutesreakdnpath} --mu_tesfakeup=${mutesfakeuppath} --mu_tesfakedn=${mutesfakednpath} --mu_jesup=${mujesuppath} --mu_jesdn=${mujesdnpath} --mu_metscalesofttermsup=${mumetscalesofttermsuppath} --mu_metscalesodttermsdn=${mumetscalesofttermsdnpath} --mu_metresosofttermsup=${mumetresosofttermsuppath} --mu_mueffup=${mumueffuppath} --mu_mueffdn=${mumueffdnpath} --mu_tau3pup=${mutau3puppath} --mu_tau3pdn=${mutau3pdnpath} --mu_mutrigup=${mumutriguppath} --mu_mutrigdn=${mumutrigdnpath} --mu_tauismuonup=${mutauismuonuppath} --mu_tauismuondn=${mutauismuondnpath} --mu_tauidup=${mutauiduppath} --mu_tauiddn=${mutauiddnpath} --mu_kfactorup=${mukfactoruppath} --mu_kfactordn=${mukfactordnpath} --efolder=${epath} --e_main=${emainpath} --e_elescaleup=${eelescaleuppath} --e_elescaledn=${eelescalednpath} --e_eleresup=${eeleresuppath} --e_tesrealup=${etesrealuppath} --e_tesrealdn=${etesreakdnpath} --e_tesfakeup=${etesfakeuppath} --e_tesfakedn=${etesfakednpath} --e_jesup=${ejesuppath} --e_jesdn=${ejesdnpath} --e_metscalesofttermsup=${emetscalesofttermsuppath} --e_metscalesodttermsdn=${emetscalesofttermsdnpath} --e_metresosofttermsup=${emetresosofttermsuppath} --e_eleeffup=${eeleeffuppath} --e_eleeffdn=${eeleeffdnpath} --e_tau3pup=${etau3puppath} --e_tau3pdn=${etau3pdnpath} --e_eletrigup=${eeletriguppath} --e_eletrigdn=${eeletrigdnpath} --e_tauiselectronup=${etauiselectronuppath} --e_tauiselectrondn=${etauiselectrondnpath} --e_tauidup=${etauiduppath} --e_tauiddn=${etauiddnpath} --e_kfactorup=${ekfactoruppath} --e_kfactordn=${ekfactordnpath}
    
fi

if [ $fake == 1 ]; then

    MUHAD_FF=`ls -d $mupath/fakefactors`
    EHAD_FF=`ls -d $epath/fakefactors`

    #calc_fakefactors.py --load=$MUHAD_FF/merged --channel=2 --Wjets --dump -e
    #FF_MUHAD_Wjets=`ls -rt | tail -n 1`
    FF_MUHAD_Wjets=`ls -d $MUHAD_FF/*tauID*Wjets* | tail -n 1`

    #calc_fakefactors.py --load=$MUHAD_FF/merged --channel=2 --multijet --dump -e
    #FF_MUHAD_multijet=`ls -rt | tail -n 1`
    FF_MUHAD_multijet=`ls -d $MUHAD_FF/*tauID*multijet* | tail -n 1`

    #calc_fakefactors.py --load=$EHAD_FF/merged --channel=1 --Wjets --dump -e
    FF_EHAD_Wjets=`ls -d $EHAD_FF/*tauID*Wjets* | tail -n 1`

    #calc_fakefactors.py --load=$EHAD_FF/merged --channel=1 --multijet --dump -e
    FF_EHAD_multijet=`ls -d $EHAD_FF/*tauID*multijet* | tail -n 1`

    paper_input_fakefactor.py --muhad ${FF_MUHAD_Wjets} --ehad ${FF_EHAD_Wjets} --muhadm ${FF_MUHAD_multijet} --ehadm ${FF_EHAD_multijet}
    
fi