#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
muon_resolution_final.py

This python script prepares the final formating of plots for the Z' muon resolution studies
"""

import optparse
import ROOT

#______________________________________________________________________________
def options():
    parser = optparse.OptionParser(description="muon_resolution")
    parser.add_option('-s', '--subset', dest='subset', action='store_true', default=False)
    parser.add_option('-g', '--grid', dest='grid', action='store_true', default=False)
    parser.add_option('-l', '--logy', dest='logy', action='store_true', default=False)
    return parser.parse_args()

#______________________________________________________________________________
def main():
    ops, args = options()

    if ops.subset:
        infile = 'muon_resolution_subset_test.DYtautausubset.hist.root'
        outfile = 'muon_res_subsub'
        all = False
    else:
        infile = 'muon_resolution_test.DYtautau.hist.root'
        outfile = 'muon_res'
        all = True

    # load file
    file_init = ROOT.TFile(infile)
    # copy histograms and profiles
    h_new_mu_pT_resolution = ROOT.TProfile()
    file_init.GetObject('h_new_mu_pT_resolution', h_new_mu_pT_resolution)
    h_new_mu_iso_pT_resolution = ROOT.TProfile()
    file_init.GetObject('h_new_mu_iso_pT_resolution', h_new_mu_iso_pT_resolution)
    h_new_mu_highpt_pT_resolution = ROOT.TProfile()
    file_init.GetObject('h_new_mu_highpt_pT_resolution', h_new_mu_highpt_pT_resolution)
    h_new_mu_highpt_iso_pT_resolution = ROOT.TProfile()
    file_init.GetObject('h_new_mu_highpt_iso_pT_resolution', h_new_mu_highpt_iso_pT_resolution)
    h_new_mu_IP_3S_veto_pT_resolution = ROOT.TProfile()
    file_init.GetObject('h_new_mu_IP_3S_veto_pT_resolution', h_new_mu_IP_3S_veto_pT_resolution)
    h_new_mu_IP_3S_veto_phi_pT_resolution = ROOT.TProfile()
    file_init.GetObject('h_new_mu_IP_3S_veto_phi_pT_resolution', h_new_mu_IP_3S_veto_phi_pT_resolution)

    #g_new_mu_pT_resolution = ROOT.TH2F()
    #file_init.GetObject('g_new_mu_pT_resolution', g_new_mu_pT_resolution)
    #g_new_mu_iso_pT_resolution = ROOT.TH2F()
    #file_init.GetObject('g_new_mu_iso_pT_resolution', g_new_mu_iso_pT_resolution)
    #g_new_mu_highpt_pT_resolution = ROOT.TH2F()
    #file_init.GetObject('g_new_mu_highpt_pT_resolution', g_new_mu_highpt_pT_resolution)
    #g_new_mu_highpt_iso_pT_resolution = ROOT.TH2F()
    #file_init.GetObject('g_new_mu_highpt_iso_pT_resolution', g_new_mu_highpt_iso_pT_resolution)
    #g_new_mu_IP_pT_resolution = ROOT.TH2F()
    #file_init.GetObject('g_new_mu_IP_pT_resolution', g_new_mu_IP_pT_resolution)
    #g_new_mu_IP_3S_veto_pT_resolution = ROOT.TH2F()
    #file_init.GetObject('g_new_mu_IP_3S_veto_pT_resolution', g_new_mu_IP_3S_veto_pT_resolution)
    #g_new_mu_IP_3S_veto_phi_pT_resolution = ROOT.TH2F()
    #file_init.GetObject('g_new_mu_IP_3S_veto_phi_pT_resolution', g_new_mu_IP_3S_veto_phi_pT_resolution)

    if all:
        h_new_mu_IP_pT_resolution = ROOT.TProfile()
        file_init.GetObject('h_new_mu_IP_pT_resolution', h_new_mu_IP_pT_resolution)
        h_new_mu_IP_3S_pT_resolution = ROOT.TProfile()
        file_init.GetObject('h_new_mu_IP_3S_pT_resolution', h_new_mu_IP_3S_pT_resolution)
        g_new_mu_IP_3S_pT_resolution = ROOT.TH2F()
        file_init.GetObject('g_new_mu_IP_3S_pT_resolution', g_new_mu_IP_3S_pT_resolution)

    profiles = [(h_new_mu_pT_resolution, 1, '2012 muon selection'),
                (h_new_mu_iso_pT_resolution, 6, '2012 muon selection+iso'),
                ]
    #histtwods = [(g_new_mu_pT_resolution, '2012 muon selection'),
    #             (g_new_mu_iso_pT_resolution, '2012 muon selection+iso'),
    #             (g_new_mu_IP_pT_resolution, '2012 muon selection+d0,z0'),
    #             ]
    if all:
        profiles.append((h_new_mu_IP_pT_resolution, 4, '2012 muon selection+d0,z0'))
        profiles.append((h_new_mu_IP_3S_pT_resolution, 30, '2012 muon selection+d0,z0+3St'))
        profiles.extend([(h_new_mu_IP_3S_veto_pT_resolution, 3, '2012 muon selection+d0,z0+3St+veto'),
                         (h_new_mu_IP_3S_veto_phi_pT_resolution, 41, '2012 muon selection+d0,z0+3St+veto+phi'),
                         (h_new_mu_highpt_pT_resolution, 28, '2012 muon selection+full highpT'),
                         (h_new_mu_highpt_iso_pT_resolution, 2, '2012 muon selection+full highpT+iso'),
                         ])
        #histtwods.extend([(g_new_mu_IP_3S_pT_resolution, '2012 muon selection+d0,z0+3St'),
        #                  (g_new_mu_IP_3S_veto_pT_resolution, '2012 muon selection+d0,z0+3St+veto'),
        #                  (g_new_mu_IP_3S_veto_phi_pT_resolution, '2012 muon selection+d0,z0+3St+veto+phi'),
        #                  (g_new_mu_highpt_pT_resolution, '2012 muon selection+full highpT'),
        #                  (g_new_mu_highpt_iso_pT_resolution, '2012 muon selection+full highpT+iso'),
        #                  ])
    else:
        profiles.extend([(h_new_mu_IP_3S_veto_pT_resolution, 3, '2012 muon selection+veto'),
                         (h_new_mu_IP_3S_veto_phi_pT_resolution, 41, '2012 muon selection+veto+phi'),
                         (h_new_mu_highpt_pT_resolution, 28, '2012 muon selection+full highpT(noveto,no3ST)'),
                         (h_new_mu_highpt_iso_pT_resolution, 2, '2012 muon selection+full highpT(noveto,no3St)+iso'),
                         ])
        #histtwods.extend([(g_new_mu_IP_3S_veto_pT_resolution, '2012 muon selection+d0,z0+veto'),
        #                  (g_new_mu_IP_3S_veto_phi_pT_resolution, '2012 muon selection+d0,z0+veto+phi'),
        #                  (g_new_mu_highpt_pT_resolution, '2012 muon selection+full highpT(no3St'),
        #                  (g_new_mu_highpt_iso_pT_resolution, '2012 muon selection+full highpT(no3St)+iso'),
        #                  ])
    
    # set style for profiles
    for hist, color, title in profiles:
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetMarkerStyle(23)
        hist.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
        hist.GetYaxis().SetTitle("Resolution / 10 GeV")
        hist.SetMaximum(3.5)
        hist.SetMinimum(-2.5)

    #draw profile
    c1 = ROOT.TCanvas("c1", "Resolution", 450, 159, 600, 500)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickY(1)
    c1.SetTicks(0,1)

    if ops.grid:
        c1.SetGridy(1)
        outfile = '%s_grid' % outfile

    if ops.logy:
        c1.SetLogy(1)
        outfile = '%s_logy' % outfile

    if all:
        prof_heading = "Resolution, all high-pT cuts"
    else:
        prof_heading = "Resolution, high-pT cuts (no veto or 3 St.)"

    h_new_mu_pT_resolution.SetTitle(prof_heading)
    h_new_mu_pT_resolution.Draw("")
    for hist, color, title in profiles:
        hist.Draw("same")
    c1.Update()

    leg = ROOT.TLegend(0.15,0.7,0.4,0.85)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.02)
    for hist, color, title in profiles:
        leg.AddEntry(hist, title, 'lep')
    leg.Draw()
    c1.Update()
    c1.SaveAs('%s_profile.png' % outfile)
    c1.Clear()

    # draw zoomed in profile
    for hist, color, title in profiles:
        hist.SetMaximum(0.4)
        hist.SetMinimum(-0.4)
        if hist == h_new_mu_pT_resolution:
            hist.SetTitle("Zoomed %s" % prof_heading)
            hist.Draw()
        else:
            hist.Draw("same")
    leg.Draw()
    c1.Update()
    c1.SaveAs('%s_profile_zoom.png' % outfile)
    c1.Clear()

    # draw extra zoomed in profile
    for hist, color, title in profiles:
        hist.SetMaximum(0.14)
        hist.SetMinimum(-0.04)
        if hist == h_new_mu_pT_resolution:
            hist.SetTitle("Extra Zoomed %s" % prof_heading)
            hist.Draw()
        else:
            hist.Draw("same")
    leg.Draw()
    c1.Update()
    c1.SaveAs('%s_profile_zoomzoom.png' % outfile)
    c1.Clear()

    #draw 2D hits
    #ROOT.gStyle.SetPalette(1)
    #index = 0
    #for hist, title in histtwods:
    #    hist.SetTitle("%s; p_{T}(#mu) [GeV]; Resolution / 10 GeV" % title)
    #    hist.Draw("colz")
    #    c1.Update()
    #    c1.SaveAs("%s_2dhist_%i.png" % (outfile, index))
    #    index += 1
    #    c1.Clear()
        
    
#_______________________________________________________________________________
if __name__ == '__main__': main()
