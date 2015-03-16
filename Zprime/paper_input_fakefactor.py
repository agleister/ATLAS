#!/usr/bin/env python
"""
paper_input_fakefactor.py

prepare inputs for lephad fake factor figure(s) in Zprime paper
"""

import optparse
import math

import ROOT, rootlogon
import metaroot

#_______________________________________________________________
def options():
    parser = optparse.OptionParser(description="options")
    parser.add_option('-e', '--ehad', dest='ehad_wjets_file', default = '')
    parser.add_option('-m', '--muhad', dest='muhad_wjets_file', default = '')
    parser.add_option( '--ehadm', dest='ehad_multijet_file', default = '')
    parser.add_option( '--muhadm', dest='muhad_multijet_file', default = '')
    return parser.parse_args()

#_______________________________________________________________
def main():

    ops, args = options()

    ew_file = ROOT.TFile(ops.ehad_wjets_file)
    mw_file = ROOT.TFile(ops.muhad_wjets_file)
    em_file = ROOT.TFile(ops.ehad_multijet_file)
    mm_file = ROOT.TFile(ops.muhad_multijet_file)

    outfile = 'lephad_fakefactor_inputs.root'

    files = [ew_file, mw_file, em_file, mm_file]

    #ehad_wjet_names = [ k.GetName() for k in ew_file.GetListOfKeys() ]
    #muhad_wjet_names = [ k.GetName() for k in mw_file.GetListOfKeys() ]
    #ehad_multijet_names = [ k.GetName() for k in em_file.GetListOfKeys() ]
    #muhad_multijet_names = [ k.GetName() for k in mm_file.GetListOfKeys() ]

    prefixes = ['ehad_wjets', 'muhad_wjets', 'ehad_multijet', 'muhad_multijet']
    #names = [ehad_wjet_names, muhad_wjet_names, ehad_multijet_names, muhad_multijet_names ]

    for prefix, filename in zip(prefixes, files):

        hist1p = filename.Get('h_tau_pt_os_tau1p')
        print hist1p
        hist3p = filename.Get('h_tau_pt_os_tau3p')
        print hist3p

        #names = [ k.GetName() for k in filename.GetListOfKeys() ]
        #out_name1p = '%s_fakefactor_os_1p' % prefix
        #out_name3p = '%s_fakefactor_os_3p' % prefix

        hist1pname = '%s_fakefactor_os_1p' % prefix
        hist3pname = '%s_fakefactor_os_3p' % prefix

        hist1p.SetName(hist1pname)
        hist3p.SetName(hist3pname)

        metaroot.file.write(hist1p, outfile)
        metaroot.file.write(hist3p, outfile)

        for hist, name in [(hist1p, hist1pname), (hist3p, hist3pname)]:
            
            #Set systematic error
            h_syst = hist.Clone('%s_syst' % name)
            nbins = h_syst.GetNbinsX()
            for bin in xrange(1, nbins+1):
                h_syst.SetBinError(bin, 0.25*h_syst.GetBinContent(bin))
                
            #Set stat error and combine with syst error
            h_error = hist.Clone('%s_error' % name)
            for bin in xrange(1, nbins+1):
                syst = h_syst.GetBinError(bin)
                stat = h_error.GetBinError(bin)
                h_error.SetBinError(bin, math.sqrt(stat*stat + syst*syst))

            metaroot.file.write(h_error, outfile)

        

#_______________________________________________________________
if __name__ == '__main__': main()
