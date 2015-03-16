#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
muon_resolution_234_final.py

This python script prepares the final formating of plots for the Z' muon resolution studies
"""

import optparse
import ROOT

#______________________________________________________________________________
def options():
    parser = optparse.OptionParser(description="muon_resolution")
    #parser.add_option('-s', '--subset', dest='subset', action='store_true', default=False)
    parser.add_option('-g', '--grid', dest='grid', action='store_true', default=False)
    parser.add_option('-l', '--logy', dest='logy', action='store_true', default=False)
    parser.add_option('-z', '--zoom', dest='zoom', action='store_true', default=False)
    parser.add_option('-s', '--selection', dest='selection', type=int, default=4)
    return parser.parse_args()

#______________________________________________________________________________
def main():
    ops, args = options()

    assert ops.selection in [2,3,4,5]

    if ops.zoom:
        zoomed = 'zoom'
        mag = 10
    else:
        zoomed = ''
        mag = 1

    # load file
    infile = 'muon_resolution_6%s_hists.DY.hist.root' % zoomed
    file_init = ROOT.TFile(infile)

    #load 2d histograms, add them to the list of input histograms
    hists_2d = []
    
    g_new_mu_pT_resolution = ROOT.TH2F()
    file_init.GetObject('g_new_mu_pT_resolution', g_new_mu_pT_resolution)
    hists_2d.append((g_new_mu_pT_resolution, 'new_mu_pT_resolution'))
    g_new_mu_3S_pT_resolution = ROOT.TH2F()
    file_init.GetObject('g_new_mu_3S_pT_resolution', g_new_mu_3S_pT_resolution)
    hists_2d.append((g_new_mu_3S_pT_resolution, 'new_mu_3S_pT_resolution'))
    g_new_mu_2S3S_pT_resolution = ROOT.TH2F()
    file_init.GetObject('g_new_mu_2S3S_pT_resolution', g_new_mu_2S3S_pT_resolution)
    hists_2d.append((g_new_mu_2S3S_pT_resolution, 'new_mu_2S3S_pT_resolution'))
    g_new_mu_3S_consistent_pT_resolution = ROOT.TH2F()
    file_init.GetObject('g_new_mu_3S_consistent_pT_resolution', g_new_mu_3S_consistent_pT_resolution)
    hists_2d.append((g_new_mu_3S_consistent_pT_resolution, 'new_mu_3S_consistent_pT_resolution'))
    g_new_mu_2S3S_consistent_pT_resolution = ROOT.TH2F()
    file_init.GetObject('g_new_mu_2S3S_consistent_pT_resolution', g_new_mu_2S3S_consistent_pT_resolution)
    hists_2d.append((g_new_mu_2S3S_consistent_pT_resolution, 'new_mu_2S3S_consistent_pT_resolution'))


    for h,t  in hists_2d:
        h.SetDirectory(0)

    file_init.Close()

    file_out = ROOT.TFile('muon_resolution_projection_6%s.root' % zoomed , 'RECREATE')

    numbins = 13

    #draw profile
    c1 = ROOT.TCanvas("c1", "Resolution", 450, 159, 600, 500)
    #ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTickY(1)
    c1.SetTicks(0,2)

    # create projection histograms (of resolution split by pT range)
    for hist2d, name in hists_2d:
        hproj = [None]*numbins
        for slice in range(numbins):
            slice += 1
            #if slice == 1 and name == 'new_mu_pT_resolution':
            hproj[slice-1] = ROOT.TH1D()
            hproj[slice-1] = hist2d.ProjectionY("h_%s_py_%i" % (name,slice), slice, slice, "de")
            hproj[slice-1].Write()
            #c1.Update()
            #print 'h_%s_%i' % (name,slice)
        #h_proj.Draw()
        #c1.Update()
        
    #file_out.Close()
    
    #create list of 1dhists grouped by pT bins
    keylist = [None]*5
    for x in range(5):
        keylist[x] = x
    namelist = ['new_mu_pT_resolution','new_mu_3S_pT_resolution','new_mu_2S3S_pT_resolution']
    namelist.extend(['new_mu_3S_consistent_pT_resolution','new_mu_2S3S_consistent_pT_resolution'])
    colorlist = [1, 6, 4, 3, 2]
    legendlist = ['2012 muon selection', '2012+3st', '2012+3or2st']
    legendlist.extend(['2012+3st+IDMScon', '2012+3or2st+IDMScon'])
    
    keynamelist = zip(keylist, namelist)
    allinfolist = zip(keylist, namelist, colorlist, legendlist)
    keysize = len(keylist)

    histos = [None]*keysize
    for x in range(keysize):
        histos[x] = [None]*numbins

    for key, name in keynamelist:
        for bin in range(numbins):
            histos[key][bin] = ROOT.TH1D()
            file_out.GetObject('h_%s_py_%i' % (name, bin+1), histos[key][bin])

    #set up arrays for mean, rms, and associated errors
    means = [None]*keysize
    for x in range(keysize):
        means[x] = [None]*numbins
    rmss = [None]*keysize
    for x in range(keysize):
        rmss[x] = [None]*numbins
    meanerrors = [None]*keysize
    for x in range(keysize):
        meanerrors[x] = [None]*numbins
    rmserrors = [None]*keysize
    for x in range(keysize):
        rmserrors[x] = [None]*numbins
    tailfracs = [None]*keysize
    for x in range(keysize):
        tailfracs[x] = [None]*numbins
    

    #define some drawing options
    ROOT.gStyle.SetOptStat(0)
    if ops.logy:
        c1.SetLogy(1)
        log = '_logy'
    else:
        log = ''

    #plot similar pT bins together
    for bin in range(numbins):
        leg = ROOT.TLegend(0.60,0.7,0.85,0.85)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        leg.SetTextSize(0.02)
        for key, name, color, label in allinfolist:
            histos[key][bin].SetLineColor(color)
            histos[key][bin].SetMarkerColor(color)
            histos[key][bin].SetMarkerStyle(23)
            if key == 0:
                if bin == numbins - 1:
                    histos[key][bin].SetTitle("resolution in p_{T} range %i GeV +" % (100*bin/mag))
                else:
                    histos[key][bin].SetTitle("resolution in p_{T} range %i-%i GeV" % (100*bin/mag, 100*bin/mag +100/mag))
                histos[key][bin].GetXaxis().SetTitle("resolution")
                histos[key][bin].Draw("")
            else:
                histos[key][bin].Draw("same")
            leg.AddEntry(histos[key][bin], label, 'lep')
            #add info for stats arrays
            means[key][bin] = histos[key][bin].GetMean()
            rmss[key][bin] = histos[key][bin].GetRMS()
            meanerrors[key][bin] = histos[key][bin].GetMeanError()
            rmserrors[key][bin] = histos[key][bin].GetRMSError()
            tailfracs[key][bin] = histos[key][bin].Integral(301,501)/histos[key][bin].GetEntries()
        leg.Draw()
        c1.Update()
        c1.SaveAs('muon_resolution_6%s_pT_%i%s.pdf' % (zoomed, bin*100/mag, log))
        c1.Clear()

    #print 'mean errors %s :' % meanerrors

    #create graphs of means and rms' and tailfractions
    meangraphs = [None]*keysize
    rmsgraphs = [None]*keysize
    tailfracgraphs = [None]*keysize
    for x in range(keysize):
        meangraphs[x] = ROOT.TGraphErrors(numbins)
        rmsgraphs[x] = ROOT.TGraphErrors(numbins)
        tailfracgraphs[x] = ROOT.TGraphErrors(numbins)
    
    meanmultigraph = ROOT.TMultiGraph()
    rmsmultigraph = ROOT.TMultiGraph()
    tailfracmultigraph = ROOT.TMultiGraph()

    #fill graphs
    for key, name, color, label in allinfolist:
        meangraphs[key].SetTitle(label)
        meangraphs[key].SetLineColor(color)
        meangraphs[key].SetMarkerColor(color)
        meangraphs[key].SetMarkerStyle(8)
        rmsgraphs[key].SetTitle(label)
        rmsgraphs[key].SetLineColor(color)
        rmsgraphs[key].SetMarkerColor(color)
        rmsgraphs[key].SetMarkerStyle(8)
        tailfracgraphs[key].SetTitle(label)
        tailfracgraphs[key].SetLineColor(color)
        tailfracgraphs[key].SetMarkerColor(color)
        tailfracgraphs[key].SetMarkerStyle(8)
        for bin in range(numbins):
            meangraphs[key].SetPoint(bin, 100*bin/mag+50/mag, means[key][bin])
            meangraphs[key].SetPointError(bin, 50/mag, meanerrors[key][bin])
            rmsgraphs[key].SetPoint(bin, 100*bin/mag+50/mag, rmss[key][bin])
            rmsgraphs[key].SetPointError(bin, 50/mag, rmserrors[key][bin])
            tailfracgraphs[key].SetPoint(bin, 100*bin/mag+50/mag, tailfracs[key][bin])
            tailfracgraphs[key].SetPointError(bin, 50/mag, 0)
        meanmultigraph.Add(meangraphs[key].Clone())
        rmsmultigraph.Add(rmsgraphs[key].Clone())
        tailfracmultigraph.Add(tailfracgraphs[key].Clone())
        
    mgs = [(meanmultigraph, 'mean', 'Resolution Mean for each pT value in each subset of cuts'),
           (rmsmultigraph, 'rms', 'Resolution RMS for each pT value in each subset of cuts'),
           (tailfracmultigraph, 'tailfrac', 'Resolution Tail Fraction for each pT value in each subset of cuts')]

    for graph, name, heading in mgs:
        graph.Draw("AP")
        graph.SetTitle(heading)
        graph.GetXaxis().SetTitle("p_{T}(#mu) [GeV]")
        if name == 'tailfrac':
            graph.GetYaxis().SetTitle("Tail Fraction / 100 GeV")
        else:
            graph.GetYaxis().SetTitle("Resolution / 100 GeV")
        if ops.logy:
            c1.BuildLegend(0.50, 0.18, 0.85, 0.45)
        else:
            c1.BuildLegend(0.15, 0.58, 0.5, 0.85)
        c1.Update()
        c1.SaveAs('muon_resolution_6%s_%s_plot%s.pdf' % (zoomed, name, log))
        c1.Clear()

#_______________________________________________________________________________
if __name__ == '__main__': main()
