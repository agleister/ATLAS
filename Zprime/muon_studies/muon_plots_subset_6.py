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
        infile = 'muon_efficiency_6_hists.DY.hist.root'
        outfile = 'muon_efficiency_6'
 
    # Load histograms from file
    file_init = ROOT.TFile(infile)
    h_true_mu_pT = ROOT.TH1F()
    file_init.GetObject('h_true_mu_pT', h_true_mu_pT)
    h_new_mu_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_pT', h_new_mu_pT)
    h_new_mu_3S_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_3S_pT', h_new_mu_3S_pT)
    h_new_mu_2S3S_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_2S3S_pT', h_new_mu_2S3S_pT)
    h_new_mu_3S_consistent_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_3S_consistent_pT', h_new_mu_3S_consistent_pT)
    h_new_mu_2S3S_consistent_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_2S3S_consistent_pT', h_new_mu_2S3S_consistent_pT)

    #divide histogram by true_mu_pt to get efficiency
    g_eff_new = ROOT.TGraphAsymmErrors()
    g_eff_new_3S = ROOT.TGraphAsymmErrors()
    g_eff_new_2S3S = ROOT.TGraphAsymmErrors()
    g_eff_new_3S_consistent = ROOT.TGraphAsymmErrors()
    g_eff_new_2S3S_consistent = ROOT.TGraphAsymmErrors()
    
    g_eff_new.Divide(h_new_mu_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_3S.Divide(h_new_mu_3S_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_2S3S.Divide(h_new_mu_2S3S_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_3S_consistent.Divide(h_new_mu_3S_consistent_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_2S3S_consistent.Divide(h_new_mu_2S3S_consistent_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')

    mg_subset = ROOT.TMultiGraph()

    # set marker/line colors to distinguish distributions
    histograms = [(g_eff_new, 1, 'new selection'),
                  (g_eff_new_3S, 4, 'new selection+3st'),
                  (g_eff_new_2S3S, 30, 'new selection+3or2st'),
                  (g_eff_new_3S_consistent, 2, 'new selection+3st+consistency'),
                  (g_eff_new_2S3S_consistent, 44, 'new selection+3or2st+consistency'), ]

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

    graphs= [ (mg_subset, 'muon_highpt_subsel_noveto3st', 'Subselections of High pT Criteria'),
              ]

    for graph, name, heading in graphs:
        graph.Draw("AP")
        graph.SetTitle(heading)
        graph.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
        graph.GetYaxis().SetTitle("Efficiency / 10 GeV")
        c1.BuildLegend(0.2, 0.13, 0.6, 0.4)
        c1.Update()
        c1.SaveAs('%s.efficiency.pdf' % outfile)
        c1.Clear()
    #for hist, name  in histograms:
    #    hist.SetMarkerStyle(8)
    #    hist.Draw("AP")
    #    c1.Update()
    #    c1.SaveAs('%s.%s.png' % (outfile, name))
    

#_______________________________________________________________________________
if __name__ == '__main__': main()
