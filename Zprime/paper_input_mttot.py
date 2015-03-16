#!/usr/bin/env python
"""
paper_input_mttot.py

prepare input for lephad mttot figures in Zprime paper
"""

import argparse

import ROOT, rootlogon
import metaroot

#_______________________________________________________________
def options():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mufolder", dest="mufolder", default="")
    parser.add_argument("--efolder", dest="efolder", default="")
    parser.add_argument("--mu_main", dest="mu_main")
    parser.add_argument("--mu_muonscaleup", dest="mu_muonscaleup")
    parser.add_argument("--mu_muonscaledn", dest="mu_muonscaledn")
    parser.add_argument("--mu_muonidup", dest="mu_muonidup")
    parser.add_argument("--mu_muonmsup", dest="mu_muonmsup")
    parser.add_argument("--mu_tesrealup", dest="mu_tesrealup")
    parser.add_argument("--mu_tesrealdn", dest="mu_tesrealdn")
    parser.add_argument("--mu_tesfakeup", dest="mu_tesfakeup")
    parser.add_argument("--mu_tesfakedn", dest="mu_tesfakedn")
    parser.add_argument("--mu_jesup", dest="mu_jesup")
    parser.add_argument("--mu_jesdn", dest="mu_jesdn")
    parser.add_argument("--mu_metscalesofttermsup", dest="mu_metscalesofttermsup")
    parser.add_argument("--mu_metscalesodttermsdn", dest="mu_metscalesofttermsdn")
    parser.add_argument("--mu_metresosofttermsup", dest="mu_metresosofttermsup")
    parser.add_argument("--mu_mueffup", dest="mu_mueffup")
    parser.add_argument("--mu_mueffdn", dest="mu_mueffdn")
    parser.add_argument("--mu_tau3pup", dest="mu_tau3pup")
    parser.add_argument("--mu_tau3pdn", dest="mu_tau3pdn")
    parser.add_argument("--mu_mutrigup", dest="mu_mutrigup")
    parser.add_argument("--mu_mutrigdn", dest="mu_mutrigdn")
    parser.add_argument("--mu_tauismuonup", dest="mu_tauismuonup")
    parser.add_argument("--mu_tauismuondn", dest="mu_tauismuondn")
    parser.add_argument("--mu_tauidup", dest="mu_tauidup")
    parser.add_argument("--mu_tauiddn", dest="mu_tauiddn")
    parser.add_argument("--mu_kfactorup", dest="mu_kfactorup")
    parser.add_argument("--mu_kfactordn", dest="mu_kfactordn")
    parser.add_argument("--e_main", dest="e_main")
    parser.add_argument("--e_elescaleup", dest="e_elescaleup")
    parser.add_argument("--e_elescaledn", dest="e_elescaledn")
    parser.add_argument("--e_eleresup", dest="e_eleresup")
    parser.add_argument("--e_tesrealup", dest="e_tesrealup")
    parser.add_argument("--e_tesrealdn", dest="e_tesrealdn")
    parser.add_argument("--e_tesfakeup", dest="e_tesfakeup")
    parser.add_argument("--e_tesfakedn", dest="e_tesfakedn")
    parser.add_argument("--e_jesup", dest="e_jesup")
    parser.add_argument("--e_jesdn", dest="e_jesdn")
    parser.add_argument("--e_metscalesofttermsup", dest="e_metscalesofttermsup")
    parser.add_argument("--e_metscalesodttermsdn", dest="e_metscalesofttermsdn")
    parser.add_argument("--e_metresosofttermsup", dest="e_metresosofttermsup")
    parser.add_argument("--e_eleeffup", dest="e_eleeffup")
    parser.add_argument("--e_eleeffdn", dest="e_eleeffdn")
    parser.add_argument("--e_tau3pup", dest="e_tau3pup")
    parser.add_argument("--e_tau3pdn", dest="e_tau3pdn")
    parser.add_argument("--e_eletrigup", dest="e_eletrigup")
    parser.add_argument("--e_eletrigdn", dest="e_eletrigdn")
    parser.add_argument("--e_tauiselectronup", dest="e_tauiselectronup")
    parser.add_argument("--e_tauiselectrondn", dest="e_tauiselectrondn")
    parser.add_argument("--e_tauidup", dest="e_tauidup")
    parser.add_argument("--e_tauiddn", dest="e_tauiddn")
    parser.add_argument("--e_kfactorup", dest="e_kfactorup")
    parser.add_argument("--e_kfactordn", dest="e_kfactordn")

    #parser.add_argument("--MZDYtautau", dest="MZDYtautau")
    #parser.add_argument("--MWjets", dest="MWjets")
    #parser.add_argument("--MZmumu", dest="MZmumu")
    #parser.add_argument("--Mtop", dest="Mtop")
    #parser.add_argument("--Mdiboson", dest="Mdiboson")
    #parser.add_argument("--Mdata", dest="Mdata")
    #parser.add_argument("--MZp1000SSM", dest="MZp1000SSM")
    #parser.add_argument("--EZtautau", dest="EZtautau")
    #parser.add_argument("--EWjets", dest="EWjets")
    #parser.add_argument("--EZee", dest="EZee")
    #parser.add_argument("--Etop", dest="Etop")
    #parser.add_argument("--Ediboson", dest="Ediboson")
    #parser.add_argument("--Edata", dest="Edata")
    #parser.add_argument("--EZp1000SSM", dest="EZp1000SSM")
    return parser.parse_args()

#_______________________________________________________________
def main():
    
    ops = vars(options())

    outfile = 'lephad_mttot_inputs.root'
    outfiles = {}
    outfiles['e']      = 'lephad_kinematic_inputs_nosfm.root'
    outfiles['e_sfm']  = 'ehad_kinematic_inputs_sfm.root'
    outfiles['mu']     = 'lephad_kinematic_inputs_nosfm.root'
    outfiles['mu_sfm'] = 'muhad_kinematic_inputs_sfm.root'

    #files = {'MZDYtautau' : ops.MZDYtautau,
    #         'MWjets' : ops.MWJets,
    #         'MZmumu' : ops.MW

    leps = []
    channels = []

    if ops['efolder'] != "" :
        leps.append('e')
        channels.append('e')
    if ops['mufolder'] != "" :
        leps.append('m')
        channels.append('mu')

    for lep, channel in zip(leps, channels):

        if lep == 'e':
            folders = ['main', 'elescaleup', 'elescaledn', 
                       'eleresup', 'tesrealup', 'tesrealdn',
                       'tesfakeup', 'tesfakedn', 'jesup', 'jesdn',
                       'metscalesofttermsup','metscalesofttermsdn', 'metresosofttermsup', 
                       'eleeffup', 'eleeffdn', 'tau3pup', 'tau3pdn',
                       'eletrigup', 'eletrigdn', 'tauiselectronup', 'tauiselectrondn',
                       'tauidup', 'tauiddn', 'kfactorup', 'kfactordn']
        elif lep == 'm':
            folders = ['main', 'muonscaleup', 'muonscaledn', 
                       'muonidup', 'muonmsup', 'tesrealup', 'tesrealdn',
                       'tesfakeup', 'tesfakedn', 'jesup', 'jesdn',
                       'metscalesofttermsup','metscalesofttermsdn', 'metresosofttermsup', 
                       'mueffup', 'mueffdn', 'tau3pup', 'tau3pdn',
                       'mutrigup', 'mutrigdn', 'tauismuonup', 'tauismuondn',
                       'tauidup', 'tauiddn', 'kfactorup', 'kfactordn']

        #folders = ['main']
        for folder in folders:

            path_to_hists = '%s/merged' % (ops['%s_%s' % (channel, folder)])
            samples = ['Ztautau', 'Z%s%s' % (channel,channel) , 'diboson', 'top']
            #samples = ['Ztautau', 'Z%s%s' % (channel,channel) , 'diboson', 'top', 
            #           'ZprimeLH500', 'ZprimeLH625', 'ZprimeLH750', 'ZprimeLH875',
            #           'ZprimeLH1000', 'ZprimeLH1125', 'ZprimeLH1250', 'ZprimeLH1375',
            #           'ZprimeLH1500', 'ZprimeLH1625', 'ZprimeLH1750', 'ZprimeLH1875',
            #           'ZprimeLH2000', 'ZprimeLH2125', 'ZprimeLH2250', 'ZprimeLH2375', 'ZprimeLH2500',
            #           'ZprimeRH500', 'ZprimeRH625', 'ZprimeRH750', 'ZprimeRH875',
            #           'ZprimeRH1000', 'ZprimeRH1125', 'ZprimeRH1250', 'ZprimeRH1375',
            #           'ZprimeRH1500', 'ZprimeRH1625', 'ZprimeRH1750', 'ZprimeRH1875',
            #           'ZprimeRH2000', 'ZprimeRH2125', 'ZprimeRH2250', 'ZprimeRH2375', 'ZprimeRH2500',
            #           'ZprimeSSM500', 'ZprimeSSM625', 'ZprimeSSM750', 'ZprimeSSM875',
            #           'ZprimeSSM1000', 'ZprimeSSM1125', 'ZprimeSSM1250', 'ZprimeSSM1375',
            #           'ZprimeSSM1500', 'ZprimeSSM1625', 'ZprimeSSM1750', 'ZprimeSSM1875',
            #           'ZprimeSSM2000', 'ZprimeSSM2125', 'ZprimeSSM2250', 'ZprimeSSM2375', 'ZprimeSSM2500']
            for signal in ['LH', 'RH', 'SSM', 'Narrow', 'Wide']:
                for i in xrange(17):
                    samples.append('Zprime%s%i' % (signal, 500 + i*125))
            for i in xrange(17):
                for s in [3, 10, 20, 30, 40, 50, 60, 70, 80, 90, 96]:
                    samples.append('ZprimeSFM%i_s2p%02d' % (500 + i*125, s))
            if folder == 'main':
                samples.extend(['WjetsFakeTaus', 'mydata'])

            #samples = ['Ztautau']
            for sample in samples:
                if "SFM" in sample:
                    sfm = "_sfm"
                else:
                    sfm = ""
                infile = ROOT.TFile('%s/%s/hadd.root' % (path_to_hists, sample))
                print infile
                if infile:
                    hist_dir_evt = 'regions/Zprime/lepPt30_tauPt30_tauLepVeto_os_tauLepDphi_mT50/event'
                    hist_dir_tau = 'regions/Zprime/lepPt30_tauPt30_tauLepVeto_os_tauLepDphi_mT50/tau'
                    hist_dir_lep = 'regions/Zprime/lepPt30_tauPt30_tauLepVeto_os_tauLepDphi_mT50/lep'
                    hist_dir_met = 'regions/Zprime/lepPt30_tauPt30_tauLepVeto_os_tauLepDphi_mT50/met'
                    hist_mttot_1  = infile.Get('%s/h_mTtot_varbins' % hist_dir_evt)
                    hist_mttot_2  = infile.Get('%s/h_mTtot_varbins2' % hist_dir_evt)
                    hist_pt_tau   = infile.Get('%s/h_tau_pt_varbins' % hist_dir_tau)
                    hist_pt_tau_2 = infile.Get('%s/h_tau_pt_varbins2' % hist_dir_tau)
                    hist_numtrack_tau = infile.Get('%s/h_tau_numTrack' % hist_dir_tau)
                    hist_pt_lep   = infile.Get('%s/h_lep_pt_varbins' % hist_dir_lep)
                    hist_pt_lep_2 = infile.Get('%s/h_lep_pt_varbins2' % hist_dir_lep)
                    hist_met      = infile.Get('%s/h_met_varbins' % hist_dir_met)
                    hist_mt       = infile.Get('regions/Zprime/lepPt30_tauPt30_tauLepVeto_os/event/h_mT')
                    hist_mt_dphi  = infile.Get('regions/Zprime/lepPt30_tauPt30_tauLepVeto_os_tauLepDphi/event/h_mT')
                    #hist3 = infile.Get('%s/h_mTtot_varbins3' % hist_dir)
                    #hist1.SetName('test')
                    hists = [hist_mttot_1, hist_mttot_2, hist_pt_tau, hist_pt_tau_2, 
                             hist_numtrack_tau, hist_pt_lep, hist_pt_lep_2, hist_met]
                    hist_names = ['mTtot_varbins', 'mTtot_varbins2', 'tau_pt_varbins', 'tau_pt_varbins2', 
                                  'tau_numTrack', 'lep_pt_varbins', 'lep_pt_varbins2', 'met_varbins']
                    if (not ('Zprime' in sample)) or ('ZprimeSSM' in sample):
                        hists.append(hist_mt)
                        hists.append(hist_mt_dphi)
                        hist_names.append('mT')
                        hist_names.append('mT_after_dphi')
                    for x in xrange(len(hists)):
                        if hists[x]:
                            #histname = '%shad_%s_%s_mTtot_varbins%i'% (channel, folder, sample, x+1)
                            histname = '%shad_%s_%s_%s'% (channel, folder, sample, hist_names[x])
                            hists[x].SetName(histname)
                            print hists[x]
                            metaroot.file.write(hists[x], outfiles['%s%s' % (channel, sfm)], '%shad_%s_%s' % (channel, folder, sample))
                    infile.Close()
            

#_______________________________________________________________
if __name__ == '__main__': main()
