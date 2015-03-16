#!/usr/bin/env python

## std
import optparse
import time

## ROOT
import ROOT, rootlogon
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1001

## my modules
import metaroot

#_____________________________________________________________________________
def options():
    parser = optparse.OptionParser(description="options")
    parser.add_option('-f', '--file', dest='file', default='')
    parser.add_option('-c', '--channel', dest='channel', type=int, default=2)
    return parser.parse_args()

#______________________________________________________________________________
def main():

    ops, args = options()
    timestamp = time.strftime('%Y-%m-%d')

    ## get hist_names
    tfile = ROOT.TFile(ops.file)
    hist_names = [ k.GetName() for k in tfile.GetListOfKeys() ]

    hist_names_1p = []
    hist_names_3p = []

    for hn in hist_names:
        #print hn
        if hn.count('1p') and hn[:2] == 'h_':
            print hn[12:14]
            hist_names_1p.append(hn)
        elif hn.count('3p') and hn[:2] == 'h_':
            hist_names_3p.append(hn)

    #print hist_names_os
    #print hist_names_ss

    if ops.channel == 1:
        outfile = 'plot.el-tau.wjets_mt_fakefactors_%s.canv.root' % timestamp
    elif ops.channel == 2:
        outfile = 'plot.mu-tau.wjets_mt_fakefactors_%s.canv.root' % timestamp

    for list, name, label in [(hist_names_1p, '1p', '1-prong'), (hist_names_3p, '3p', '3-prong')]:

        #get histograms
        hists = []
        hist_ratios = []
        for hname in list:
            h2 = tfile.Get(hname)
            hists.append(h2)
            #if hname.count('d0GT02'):
            #    h_denom = h2.Clone("h_denom")
            #h2_ratio = h2.Clone("%s_ratio" % hname)
            #hist_ratios.append(h2_ratio)

        print hists
        
        #print hist_ratios

        #find ratios
        #for h2r in hist_ratios:
        #    if h_denom:
        #        h2r.Divide(h_denom)
        #    else:
        #        h2r.Divide(h2r)

        #name options
        plot_name = 'h_W+jet_mT_fakefactors_inclusive_%s' % name
        #ratio_name = 'ff / d0GT02'
        

        colors = [ 
            metaroot.style.black,
            metaroot.style.red,
            metaroot.style.blue,
            metaroot.style.green,
            metaroot.style.violet,
            metaroot.style.gray,
            ]

        pops = [ 
                metaroot.hist.PlotOptions(
                    marker_color = colors[i],
                    marker_size  = 1.5,
                    line_color   = colors[i], 
                    )   
                for i in xrange(len(list))
            ]   
        canvas_options = metaroot.hist.CanvasOptions()
        canvas_options_2 = metaroot.hist.CanvasOptions(width=900, height=900)

        draw_options = ['PE']*len(hists)

        title = ';%s;%s fakefactors' % (hists[0].GetXaxis().GetTitle(), name)
        #ratio_title = ';%s;%s' % (hists[0].GetXaxis().GetTitle(), ratio_name)

        y_min = 0.0
        y_max = 0.3

        #bin_fractional_errors = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25]
        h_syst = hists[0].Clone('h_syst')
        n_bins = h_syst.GetNbinsX()
        #assert nbins == len(bin_fractional_errors)
        for bin in xrange(1, n_bins+1):
            #h_syst.SetBinError(bin, bin_fractional_errors[bin-1]*h_syst.GebBinContent(bin))
            h_syst.SetBinError(bin, 0.30*h_syst.GetBinContent(bin))
        hists.append(h_syst)
        #labels.append('systematic')
        draw_options.append('E2')
        pops.append( metaroot.hist.PlotOptions(fill_color=metaroot.style.yellow) )

        ## make plot
        plot = metaroot.hist.pile_hists(
            hists = hists,
            name = plot_name,
            title = title,
            draw_options = draw_options,
            plot_options = pops,
            canvas_options = canvas_options,
            min = y_min,
            max = y_max,
            )

        labels = []
        for ln in list:
            labels.append( get_label(ln) )
        labels.append('systematic')
   
        ## make legend
        legend = metaroot.hist.make_legend(
            hists = hists,
            labels = labels,
            draw_options = ['PE']*(len(hists)-1) + ['HIST'],
            x1 = 0.75,
            y1 = 0.78,
            height = 0.03,
            )
        legend.SetTextSize(0.04)
        legend.Draw()

        #make ratio plot
        
        #top_plot = plot
        #canv_top = top_plot['canvas']

        #bottom_plot = metaroot.hist.pile_hists(
        #    hists = hist_ratios,
        #    name = plot_name+'_ratios',
        #    title = ratio_title,
        #    draw_options = ['PE']*len(hists),
        #    plot_options = pops,
        #    canvas_options = canvas_options,
        #    min = 0.0,
        #    max = 2.0
        #    )
        
        #canv_bottom = bottom_plot['canvas']
        #line10 = metaroot.utils.draw_horiz_line(canv_bottom, 1.0,color=ROOT.kGray+2, width=1, style=7)

        #shared_plot = metaroot.plot.plot_shared_axis(
        #    canv_top,
        #    canv_bottom,
        #    name = plot_name+'_ratio',
        #    split = 0.35,
        #    axissep = 0.02,
        #    ndivs = [505,503],
        #    canvas_options = canvas_options_2
        #    )

        #plot = shared_plot

        print plot['canvas']
        print 'writing %s' % plot['canvas'].GetName()
        metaroot.file.write(plot['canvas'], outfile)

    print 'Done!'

#_____________________________________________________________________________
def get_label(ln):
    
    if ln.count('os'):
        label = 'os'
    else:
        label = 'ss'
    if ln.count('data'):
        label = 'data %s' % label
    elif ln.count('mc'):
        label = 'mc %s' % label
    return label

#_____________________________________________________________________________
def startswith_strip(s, x):
    if s.startswith(x):
        return s[len(x):]
    else:
        return s

#_____________________________________________________________________________
def endswith_strip(s, x):
    if s.endswith(x):
        return s[:-1*len(x)]
    else:
        return s

#_____________________________________________________________________________
if __name__ == '__main__': main()
