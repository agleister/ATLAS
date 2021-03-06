#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
muon_resolution.py

This script runs the procedure for making the plots plots used for muon resolution studies
"""

## ROOT
import ROOT
ROOT.gROOT.SetBatch(True)

## pyframe
import pyframe

## zprime 
import ZprimeTauTauLepHad2012 as Zp12

## testscripts
import TestScripts as TS

GeV = 1000.0

def main():
    version = 'DYtautausubsubset'
    input_file_dir_1 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.158736.Pythia8_AU2CTEQ6L1_DYtautau_1000M1250.merge.NTUP_TAU.e1518_s1499_s1504_r3658_r3549_p1344_tid01109459_00'
    input_file_dir_2 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.158740.Pythia8_AU2CTEQ6L1_DYtautau_2000M2250.merge.NTUP_TAU.e1518_s1499_s1504_r3658_r3549_p1344_tid01109463_00'
    input_file_dir_3 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.158738.Pythia8_AU2CTEQ6L1_DYtautau_1500M1750.merge.NTUP_TAU.e1518_s1499_s1504_r3658_r3549_p1344_tid01109461_00'
    input_file_dir_4 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.158741.Pythia8_AU2CTEQ6L1_DYtautau_2250M2500.merge.NTUP_TAU.e1518_s1499_s1504_r3658_r3549_p1344_tid01109464_00'
    
    input_files = [
        '%s/NTUP_TAU.01109459._000001.root.1' % input_file_dir_1,
        '%s/NTUP_TAU.01109459._000002.root.1' % input_file_dir_1,
        '%s/NTUP_TAU.01109459._000003.root.1' % input_file_dir_1,
        '%s/NTUP_TAU.01109459._000004.root.1' % input_file_dir_1,
        '%s/NTUP_TAU.01109459._000005.root.2' % input_file_dir_1,
        '%s/NTUP_TAU.01109459._000006.root.2' % input_file_dir_1,
        '%s/NTUP_TAU.01109463._000001.root.1' % input_file_dir_2,
        #'%s/NTUP_TAU.01109463._000002.root.1' % input_file_dir_2,
        '%s/NTUP_TAU.01109463._000003.root.1' % input_file_dir_2,
        '%s/NTUP_TAU.01109463._000004.root.1' % input_file_dir_2,
        '%s/NTUP_TAU.01109463._000005.root.1' % input_file_dir_2,
        '%s/NTUP_TAU.01109461._000001.root.1' % input_file_dir_3,
        '%s/NTUP_TAU.01109461._000002.root.1' % input_file_dir_3,
        '%s/NTUP_TAU.01109461._000003.root.1' % input_file_dir_3,
        '%s/NTUP_TAU.01109461._000004.root.1' % input_file_dir_3,
        '%s/NTUP_TAU.01109461._000005.root.1' % input_file_dir_3,
        '%s/NTUP_TAU.01109461._000006.root.1' % input_file_dir_3,
        '%s/NTUP_TAU.01109464._000001.root.1' % input_file_dir_4,
        '%s/NTUP_TAU.01109464._000002.root.1' % input_file_dir_4,
        '%s/NTUP_TAU.01109464._000003.root.1' % input_file_dir_4,
        '%s/NTUP_TAU.01109464._000004.root.1' % input_file_dir_4,
        '%s/NTUP_TAU.01109464._000005.root.1' % input_file_dir_4,
        '%s/NTUP_TAU.01109464._000006.root.1' % input_file_dir_4,
        ]

        #build the chain_________________________________
    chain = ROOT.TChain('tau')
    for file in input_files:
        chain.Add(file)

    #configure event loop
    loop = pyframe.core.EventLoop(name='muon_resolution_test', version=version, quiet=False)

    #build vertices
    Zp12.algs.BuildAlg.BuildVertices(loop)

    #build muon objects
    prefixes = []
    keys = []
    
    prefixes.append('mu_staco_') # staco muons
    keys.append('staco_muons')

    loop += pyframe.algs.ListBuilder(prefixes = prefixes, keys = keys)
    loop += pyframe.p4calc.AttachTLVs(keys = keys)

    #corrections
    loop += pyframe.muon.MuonMomentumCorrections(key='staco_muons')
    loop += pyframe.muon.MuonIsolationVariables(   key='staco_muons')
    loop += pyframe.muon.MuonIsolationCorrections( key='staco_muons')

    #define muon selections
    #new muon selection (without high pT)
    muon_new_selector = pyframe.muon.MuonSelector(
        min_pt = 4.0*GeV,
        max_abs_eta = 2.5,
        flags = ['loose', 'isCombinedMuon', 'truth_matched'],
        min_expected_pixel_hits = 1,
        min_expected_sct_hits = 5,
        max_si_holes = 2,
        req_trt_quality_new = True,
        )
    loop += pyframe.selectors.SelectorAlg('NewMuonSelector',
                                          selector = muon_new_selector,
                                          key_in='staco_muons',
                                          key_out='new_muons'
                                          )
    muon_new_iso_selector = pyframe.muon.MuonSelector(
        max_scores = [('ptcone30rel', 0.05),
                      ],
        )
    #add criteria for high pT one-by-one
    muon_new_IP_selector = pyframe.muon.MuonSelector(
        max_d0 = 0.2,
        max_z0 = 1.0,   
        )
    #muon_new_3S_selector = pyframe.muon.MuonSelector(
    #    highpt_three_station = True,
    #    )
    #muon_new_veto_selector = pyframe.muon.MuonSelector(
    #    highpt_veto_MS_chambers = True,
    #    )
    muon_new_phi_selector = pyframe.muon.MuonSelector(
        highpt_phi_hits = True,
        )
    muon_new_IDMS_selector = pyframe.muon.MuonSelector(
        ID_MS_consistency = True,
        )
    #create lists of muons based on criteria selection
    loop += pyframe.selectors.SelectorAlg('NewMuonNewIsoSelector',
                                          selector = muon_new_iso_selector,
                                          key_in='new_muons',
                                          key_out='new_iso_muons',
                                          )
    loop += pyframe.selectors.SelectorAlg('NewMuonIPSelector',
                                        selector = muon_new_IP_selector,
                                        key_in='new_muons',
                                        key_out='new_IP_muons',
                                        )
    #loop += pyframe.selectors.SelectorAlg('NewMuon3sSelector',
    #                                    selector = muon_new_3S_selector,
    #                                    key_in='new_IP_muons',
    #                                    key_out='new_IP_3S_muons',
    #                                      )
    #loop += pyframe.selectors.SelectorAlg('NewMuonVetoSelector',
    #                                    selector = muon_new_veto_selector,
    #                                    key_in='new_IP_muons',
    #                                    key_out='new_IP_3S_veto_muons',
    #                                    )
    loop += pyframe.selectors.SelectorAlg('NewMuonPhiSelector',
                                        selector = muon_new_phi_selector,
                                        key_in='new_IP_muons',
                                        key_out='new_IP_3S_veto_phi_muons',
                                        )
    loop += pyframe.selectors.SelectorAlg('NewMuonIDMSSelector',
                                        selector = muon_new_IDMS_selector,
                                        key_in='new_IP_3S_veto_phi_muons',
                                        key_out='new_highpt_muons',
                                        )
    loop += pyframe.selectors.SelectorAlg('NewMuonHighPtIsoSelector',
                                        selector = muon_new_iso_selector,
                                        key_in='new_highpt_muons',
                                        key_out='new_highpt_iso_muons',
                                        )
    #make plots
    loop += pyframe.algs.LooperAlg(
        key = 'new_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_',
        dir='',
        )
    loop += pyframe.algs.LooperAlg(
        key = 'new_iso_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_iso_',
        dir='',
        )
    loop += pyframe.algs.LooperAlg(
        key = 'new_IP_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_IP_',
        dir='',
        )
    #loop += pyframe.algs.LooperAlg(
    #    key = 'new_IP_3S_muons',
    #    func = TS.muon_resolution_plots.plot_mu_pt,
    #    prefix = 'new_mu_IP_3S_',
    #    dir='',
    #    )
    #loop += pyframe.algs.LooperAlg(
    #    key = 'new_IP_3S_veto_muons',
    #    func = TS.muon_resolution_plots.plot_mu_pt,
    #    prefix = 'new_mu_IP_3S_veto_',
    #    dir='',
    #    )
    loop += pyframe.algs.LooperAlg(
        key = 'new_IP_3S_veto_phi_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_IP_3S_veto_phi_',
        dir='',
        )
    loop += pyframe.algs.LooperAlg(
        key = 'new_highpt_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_highpt_',
        dir='',
        )
    loop += pyframe.algs.LooperAlg(
        key = 'new_highpt_iso_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_highpt_iso_',
        dir='',
        )

    #run the chain
    loop.run(chain, 0, -1)

#__________________________________________________________________
if __name__ == '__main__': main()
    
    

                                     
