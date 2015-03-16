#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
muon_plot_final.py

This program makes the final efficiency plots for the muon selection studies
"""

import optparse
import ROOT

#______________________________________________________________________________
def options():
    parser = optparse.OptionParser(description="muon_efficiency")
    parser.add_option('-d', '--data', dest='data', action='store_true', default=False)
    return parser.parse_args()

#______________________________________________________________________________
def main():

    ops, args = options()

    if ops.data:
        infile = 'muon_efficiency_test.data_period_A.hist.root'
        outfile = 'data_period_A'
    else:
        infile = 'muon_efficiency_test.DYtautau1000M1250.hist.root'
        outfile = 'DYtautau'

    # Load histograms from file
    file_init = ROOT.TFile(infile)
    h_new_mu_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_pT', h_new_mu_pT)
    h_true_mu_pT = ROOT.TH1F()
    file_init.GetObject('h_true_mu_pT', h_true_mu_pT)
    h_old_mu_pT = ROOT.TH1F()
    file_init.GetObject('h_old_mu_pT', h_old_mu_pT)
    h_new_highpt_mu_pT = ROOT.TH1F()
    file_init.GetObject('h_new_highpt_mu_pT', h_new_highpt_mu_pT)
    h_old_mu_old_iso_pT = ROOT.TH1F()
    file_init.GetObject('h_old_mu_old_iso_pT', h_old_mu_old_iso_pT)
    h_new_mu_old_iso_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_old_iso_pT', h_new_mu_old_iso_pT)
    h_new_mu_new_iso_pT = ROOT.TH1F()
    file_init.GetObject('h_new_mu_new_iso_pT', h_new_mu_new_iso_pT)
    h_new_highpt_mu_iso_pT = ROOT.TH1F()
    file_init.GetObject('h_new_highpt_mu_iso_pT', h_new_highpt_mu_iso_pT)

    #divide histogram by true_mu_pt to get efficiency
    g_eff_old = ROOT.TGraphAsymmErrors()
    g_eff_new = ROOT.TGraphAsymmErrors()
    g_eff_new_highpt = ROOT.TGraphAsymmErrors()
    g_eff_old_oldiso = ROOT.TGraphAsymmErrors()
    g_eff_new_oldiso = ROOT.TGraphAsymmErrors()
    g_eff_new_newiso = ROOT.TGraphAsymmErrors()
    g_eff_new_highptiso = ROOT.TGraphAsymmErrors()
    g_eff_old.Divide(h_old_mu_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new.Divide(h_new_mu_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_highpt.Divide(h_new_highpt_mu_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_old_oldiso.Divide(h_old_mu_old_iso_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_oldiso.Divide(h_new_mu_old_iso_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_newiso.Divide(h_new_mu_new_iso_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')
    g_eff_new_highptiso.Divide(h_new_highpt_mu_iso_pT,h_true_mu_pT,'cl=0.683 b(1,1) mode')

    # set marker/line colors to distinguish distributions
    histograms = [(g_eff_old, 1, 'old selection, no isolation'),
                  (g_eff_new, 4, 'new selection, no isolation'),
                  (g_eff_new_highpt, 2, 'new high-pT selection, no isolation'),
                  (g_eff_old_oldiso, 15, 'old selection, old isolation'),
                  (g_eff_new_oldiso, 7, 'new selection, old isolation'),
                  (g_eff_new_newiso, 38, 'new selection, new isolation'),
                  (g_eff_new_highptiso, 28, 'new high-Pt selection, new isolation'),]

    for hist, color, title  in histograms:
        hist.SetTitle(title)
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetMarkerStyle(8)
        hist.GetYaxis().SetLimits(0., 1.05)
        #hist.GetHistogram().SetMaximum(1.05)

    # set up multigraphs
    
    mg_old_all = ROOT.TMultiGraph()
    mg_new_all = ROOT.TMultiGraph()
    mg_new_highpt_all = ROOT.TMultiGraph()
    mg_no_iso = ROOT.TMultiGraph()
    mg_all_iso = ROOT.TMultiGraph()
    mg_all = ROOT.TMultiGraph()

    mg_old_all.Add(g_eff_old.Clone())
    mg_old_all.Add(g_eff_old_oldiso.Clone())
    
    
    mg_new_all.Add(g_eff_new.Clone())
    mg_new_all.Add(g_eff_new_oldiso.Clone())
    mg_new_all.Add(g_eff_new_newiso.Clone())

    mg_new_highpt_all.Add(g_eff_new_highpt.Clone())
    mg_new_highpt_all.Add(g_eff_new_highptiso.Clone())
    
    mg_no_iso.Add(g_eff_old.Clone())
    mg_no_iso.Add(g_eff_new.Clone())
    mg_no_iso.Add(g_eff_new_highpt.Clone())
    
    mg_all_iso.Add(g_eff_old_oldiso.Clone())
    mg_all_iso.Add(g_eff_new_oldiso.Clone())
    mg_all_iso.Add(g_eff_new_newiso.Clone())
    mg_all_iso.Add(g_eff_new_highptiso.Clone())

    mg_all.Add(g_eff_old.Clone())
    mg_all.Add(g_eff_new.Clone())
    mg_all.Add(g_eff_new_highpt.Clone())
    mg_all.Add(g_eff_old_oldiso.Clone())
    mg_all.Add(g_eff_new_oldiso.Clone())
    mg_all.Add(g_eff_new_newiso.Clone())
    mg_all.Add(g_eff_new_highptiso.Clone())
    
    #Set up legends for graphs
    #leg_old = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
    #leg_new = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
    #leg_new_highpt = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
    #leg_no_iso = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
    #leg_all_iso = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)
    #leg_all = ROOT.TLegend(0.5, 0.7, 0.9, 0.9)

    #legends = [(leg_old, [g_effold]
    #leg_old.Set
    #leg_old.SetHeader("Efficiency")

    #histograms = [(g_eff_old,'h_old_mu_pT'),
    #              (g_eff_new, 'h_new_mu_pT'),
    #              (g_eff_new_highpt, 'h_new_highpt_mu_pT'),
    #              (g_eff_old_oldiso, 'h_old_mu_old_iso_pT'),
    #              (g_eff_new_oldiso, 'h_new_mu_old_iso_pT'),
    #              (g_eff_new_newiso, 'h_new_mu_new_iso_pT'),
    #              (g_eff_new_highptiso, 'h_new_highpt_mu_iso_pT'),]

    c1 = ROOT.TCanvas("c1", "Efficiency of cuts", 450, 159, 600, 500)
    c1.SetGrid()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickY(1)

    #mg_old_all.Draw("AP")
    #g_eff_old.Draw("AP")
    #c1.BuildLegend(0.2, 0.13, 0.6, 0.4)
    #c1.Update()
    #c1.SaveAs("ThisIsATest.png")

    graphs= [ (mg_old_all, 'muon_old_plot', 'Old Muon Selection'),
              (mg_new_all, 'muon_new_plot', 'New Muon Selection (not high pT)'),
              (mg_new_highpt_all, 'muon_new_highpt_plot', 'New Muon Selection With High pT'),
              (mg_no_iso, 'no_isolation_plot', 'No Isolation Applied'),
              (mg_all_iso, 'all_isolation_plot', 'Isolation Applied'),
              (mg_all, 'all_plot', 'All Muon Selections'),
              ]

    for graph, name, heading in graphs:
        graph.Draw("AP")
        graph.SetTitle(heading)
        graph.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
        graph.GetYaxis().SetTitle("Efficiency / 10 GeV")
        c1.BuildLegend(0.2, 0.13, 0.6, 0.4)
        c1.Update()
        c1.SaveAs('%s.efficiency_%s.png' % (outfile, name))
        c1.Clear()
    #for hist, name  in histograms:
    #    hist.SetMarkerStyle(8)
    #    hist.Draw("AP")
    #    c1.Update()
    #    c1.SaveAs('%s.%s.png' % (outfile, name))
    

#_______________________________________________________________________________
if __name__ == '__main__': main()
