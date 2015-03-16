#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
muon_resolution6.py

This script runs the procedure for making the plots used for muon resolution studies (final cuts from Will proposal)
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
    version = 'DY'

    #Dytautau files
    input_file_dir_1 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.158736.Pythia8_AU2CTEQ6L1_DYtautau_1000M1250.merge.NTUP_TAU.e1518_s1499_s1504_r3658_r3549_p1344_tid01109459_00'
    input_file_dir_2 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.158740.Pythia8_AU2CTEQ6L1_DYtautau_2000M2250.merge.NTUP_TAU.e1518_s1499_s1504_r3658_r3549_p1344_tid01109463_00'
    input_file_dir_3 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.158738.Pythia8_AU2CTEQ6L1_DYtautau_1500M1750.merge.NTUP_TAU.e1518_s1499_s1504_r3658_r3549_p1344_tid01109461_00'
    input_file_dir_4 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.158741.Pythia8_AU2CTEQ6L1_DYtautau_2250M2500.merge.NTUP_TAU.e1518_s1499_s1504_r3658_r3549_p1344_tid01109464_00'
    #DYmumu files
    input_file_dir_5 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.129530.PowhegPythia8_AU2CT10_DYmumu_1000M1250.merge.NTUP_TAU.e1248_s1469_s1470_r3752_r3549_p1344_tid01109445_00'
    input_file_dir_6 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.129532.PowhegPythia8_AU2CT10_DYmumu_1500M1750.merge.NTUP_TAU.e1248_s1469_s1470_r3542_r3549_p1344_tid01109447_00'
    input_file_dir_7 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.129534.PowhegPythia8_AU2CT10_DYmumu_2000M2250.merge.NTUP_TAU.e1248_s1469_s1470_r3542_r3549_p1344_tid01109449_00'
    input_file_dir_8 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.129536.PowhegPythia8_AU2CT10_DYmumu_2500M2750.merge.NTUP_TAU.e1248_s1469_s1470_r3542_r3549_p1344_tid01109451_00'
    input_file_dir_9 = '/group/atlas/data/tauD3PD/zprimetautau/D3PD/mc12_8TeV.129537.PowhegPythia8_AU2CT10_DYmumu_2750M3000.merge.NTUP_TAU.e1248_s1469_s1470_r3542_r3549_p1344_tid01109452_00'
    
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
        '%s/NTUP_TAU.01109445._000001.root.1' % input_file_dir_5,
        '%s/NTUP_TAU.01109445._000002.root.1' % input_file_dir_5,
        '%s/NTUP_TAU.01109445._000003.root.1' % input_file_dir_5,
        '%s/NTUP_TAU.01109445._000004.root.1' % input_file_dir_5,
        '%s/NTUP_TAU.01109447._000001.root.1' % input_file_dir_6,
        '%s/NTUP_TAU.01109447._000002.root.1' % input_file_dir_6,
        '%s/NTUP_TAU.01109447._000003.root.1' % input_file_dir_6,
        '%s/NTUP_TAU.01109447._000004.root.1' % input_file_dir_6,
        '%s/NTUP_TAU.01109449._000001.root.1' % input_file_dir_7,
        '%s/NTUP_TAU.01109449._000002.root.1' % input_file_dir_7,
        '%s/NTUP_TAU.01109449._000003.root.1' % input_file_dir_7,
        '%s/NTUP_TAU.01109449._000004.root.1' % input_file_dir_7,
        '%s/NTUP_TAU.01109451._000001.root.1' % input_file_dir_8,
        '%s/NTUP_TAU.01109451._000002.root.4' % input_file_dir_8,
        '%s/NTUP_TAU.01109451._000003.root.2' % input_file_dir_8,
        '%s/NTUP_TAU.01109451._000004.root.2' % input_file_dir_8,
        '%s/NTUP_TAU.01109452._000001.root.1' % input_file_dir_9,
        '%s/NTUP_TAU.01109452._000002.root.1' % input_file_dir_9,
        '%s/NTUP_TAU.01109452._000003.root.1' % input_file_dir_9,
        '%s/NTUP_TAU.01109452._000004.root.1' % input_file_dir_9,
        ]

        #build the chain_________________________________
    chain = ROOT.TChain('tau')
    for file in input_files:
        chain.Add(file)

    #configure event loop
    loop = pyframe.core.EventLoop(name='muon_resolution_6_hists', version=version, quiet=False)

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
    
    #add criteria for high pT one-by-one
    
    muon_new_3S_selector = pyframe.muon.MuonSelector(
        highpt_three_station = True,
        )
    muon_new_2S_selector = pyframe.muon.MuonSelector(
        anti_highpt_three_station = True,
        highpt_two_station = True,
        tower_veto_2s = True,
        )
    muon_new_IDMS_selector = pyframe.muon.MuonSelector(
        ID_MS_consistency = True,
        )
    muon_new_IDMS_selector_2s = pyframe.muon.MuonSelector(
        ID_MS_consistency2 = True,
        )
    #create lists of muons based on criteria selection
    loop += pyframe.selectors.SelectorAlg('NewMuon3sSelector',
                                        selector = muon_new_3S_selector,
                                        key_in='new_muons',
                                        key_out='new_3S_muons',
                                          )
    loop += pyframe.selectors.SelectorAlg('NewMuon2sSelector',
                                        selector = muon_new_2S_selector,
                                        key_in='new_muons',
                                        key_out='new_2S_muons',
                                          )
    loop += pyframe.algs.ListMerger( keys_in = ['new_3S_muons', 'new_2S_muons'],
                                     key_out = 'new_2S3S_muons',
                                     )
    loop += pyframe.selectors.SelectorAlg('NewMuonIDMSSelector',
                                        selector = muon_new_IDMS_selector,
                                        key_in='new_3S_muons',
                                        key_out='new_3S_consistent_muons',
                                        )
    loop += pyframe.selectors.SelectorAlg('NewMuonIDMSSelector2s',
                                        selector = muon_new_IDMS_selector_2s,
                                        key_in='new_2S_muons',
                                        key_out='new_2S_consistent_muons',
                                        )
    loop += pyframe.algs.ListMerger( keys_in = ['new_3S_consistent_muons', 'new_2S_consistent_muons'],
                                     key_out = 'new_2S3S_consistent_muons',
                                     )
    
    #make plots
    loop += pyframe.algs.LooperAlg(
        key = 'new_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_',
        dir='',
        )
    loop += pyframe.algs.LooperAlg(
        key = 'new_3S_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_3S_',
        dir='',
        )
    loop += pyframe.algs.LooperAlg(
        key = 'new_2S3S_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_2S3S_',
        dir='',
        )
    loop += pyframe.algs.LooperAlg(
        key = 'new_3S_consistent_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_3S_consistent_',
        dir='',
        )
    loop += pyframe.algs.LooperAlg(
        key = 'new_2S3S_consistent_muons',
        func = TS.muon_resolution_plots.plot_mu_pt,
        prefix = 'new_mu_2S3S_consistent_',
        dir='',
        )

    #run the chain
    loop.run(chain, 0, -1)

#__________________________________________________________________
if __name__ == '__main__': main()
    
