#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
muon_plot_subset.py

This program makes the final efficiency plots for the muon selection studies
"""

import optparse
import ROOT

#______________________________________________________________________________
def options():
    parser = optparse.OptionParser(description="muon_efficiency")
    parser.add_option('-d', '--data', dest='data', action='store_true', default=False)
    parser.add_option('-s', '--selection', dest='selection', type=int, default=4)
    return parser.parse_args()

#______________________________________________________________________________
def main():

    ops, args = options()

    if ops.data:
        infile = 'muon_efficiency_test.data_period_A.hist.root'
        outfile = 'data_period_A'
    else:
        infile = 'muon_efficiency_%i_3S2S_test.DY.hist.root' % ops.selection
        outfile = 'muon_efficiency_%i' % ops.selection
 
    # Load histograms from file
    file_init = ROOT.TFile(infile)
    h_true_mu_pT = ROOT.TH1F()
    file_init.GetObject('h_true_mu_pT', h_true_mu_pT)
    h_new_mu_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_pT', h_new_mu_pT)
    h_new_highpt_mu_pT = ROOT.TH1F()
    file_init.GetObject('h_new_highpt_mu_pT', h_new_highpt_mu_pT)
    h_new_highpt_mu_iso_pT = ROOT.TH1F()
    file_init.GetObject('h_new_highpt_mu_iso_pT', h_new_highpt_mu_iso_pT)
    #h_new_mu_IP_pT = ROOT.TH1F()
    #file_init.GetObject('h_new_mu_IP_pT', h_new_mu_IP_pT)
    h_new_mu_IP_2S3S_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_IP_2S3S_pT', h_new_mu_IP_2S3S_pT)
    h_new_mu_IP_2S3S_veto_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_IP_2S3S_veto_pT', h_new_mu_IP_2S3S_veto_pT)
    h_new_mu_IP_2S3S_veto_phi_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_IP_2S3S_veto_phi_pT', h_new_mu_IP_2S3S_veto_phi_pT)

    #divide histogram by true_mu_pt to get efficiency
    g_eff_new = ROOT.TGraphAsymmErrors()
    #g_eff_new_IP = ROOT.TGraphAsymmErrors()
    g_eff_new_IP_2S3S = ROOT.TGraphAsymmErrors()
    g_eff_new_IP_2S3S_veto = ROOT.TGraphAsymmErrors()
    g_eff_new_IP_2S3S_veto_phi = ROOT.TGraphAsymmErrors()
    g_eff_new_highpt = ROOT.TGraphAsymmErrors()
    g_eff_new_highptiso = ROOT.TGraphAsymmErrors()
    
    g_eff_new.Divide(h_new_mu_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    #g_eff_new_IP.Divide(h_new_mu_IP_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_IP_2S3S.Divide(h_new_mu_IP_2S3S_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_IP_2S3S_veto.Divide(h_new_mu_IP_2S3S_veto_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_IP_2S3S_veto_phi.Divide(h_new_mu_IP_2S3S_veto_phi_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_highpt.Divide(h_new_highpt_mu_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_highptiso.Divide(h_new_highpt_mu_iso_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')

    mg_subset = ROOT.TMultiGraph()

    # set marker/line colors to distinguish distributions
    histograms = [(g_eff_new, 1, 'new selection'),
                  #(g_eff_new_IP, 4, 'new selection+d0,z0'),
                  (g_eff_new_IP_2S3S, 30, 'new selection+3or2st'),
                  (g_eff_new_IP_2S3S_veto, 8, 'new selection+3or2st+MCveto'),
                  (g_eff_new_IP_2S3S_veto_phi, 41, 'new selection+3or2st+MCveto+phihit'),
                  (g_eff_new_highpt, 44, 'new high-pT selection, no isolation'),
                  (g_eff_new_highptiso, 2, 'new high-Pt selection, new isolation'),]

    for hist, color, title  in histograms:
        hist.SetTitle(title)
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetMarkerStyle(8)
        hist.GetYaxis().SetLimits(0., 1.05)
        hist.GetXaxis().SetLimits(0., 1000.0)
        mg_subset.Add(hist.Clone())

    #canvas options
    c1 = ROOT.TCanvas("c1", "Efficiency of cuts", 450, 159, 600, 500)
    c1.SetGrid()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickY(1)

    #mg_old_all.Draw("AP")
    #g_eff_old.Draw("AP")
    #c1.BuildLegend(0.2, 0.13, 0.6, 0.4)
    #c1.Update()
    #c1.SaveAs("ThisIsATest.png")

    if ops.selection == 4:
        finaltitle = 'all cuts'
    elif ops.selection == 3:
        finaltitle = 'no 3 station'
    elif ops.selection == 5:
        finaltitle = 'no d0,z0'
    else:
        finaltitle = 'no IP or 3 station'

    graphs= [ (mg_subset, 'muon_highpt_subsel_noveto3st', 'Subselections of High pT Criteria (%s)' % finaltitle),
              ]

    for graph, name, heading in graphs:
        graph.Draw("AP")
        graph.SetTitle(heading)
        graph.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
        graph.GetYaxis().SetTitle("Efficiency / 10 GeV")
        c1.BuildLegend(0.2, 0.13, 0.6, 0.4)
        c1.Update()
        c1.SaveAs('%s_2S3S.efficiency.pdf' % outfile)
        c1.Clear()
    #for hist, name  in histograms:
    #    hist.SetMarkerStyle(8)
    #    hist.Draw("AP")
    #    c1.Update()
    #    c1.SaveAs('%s.%s.png' % (outfile, name))
    

#_______________________________________________________________________________
if __name__ == '__main__': main()
