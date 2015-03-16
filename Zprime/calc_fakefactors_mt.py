#!/usr/bin/env python
"""
calc_fakefactors_mt.py
"""
#------------------------------------------------------------------------------

## std
import optparse
import time

## ROOT
import ROOT, rootlogon
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1001

## my modules
import metaroot
from redict import redict

## local modules

## samples
from samples.p1344 import mc
from samples.p1443 import data

quiet = True

#_____________________________________________________________________________
def options():
    parser = optparse.OptionParser(description="options")
    parser.add_option('-d', '--dump', dest='dump', action='store_true', default=False)
    parser.add_option('-L', '--load', dest='load', type=str, default='nominal')
    parser.add_option('-c', '--channel',  dest='channel',  type=int, default=2)
    parser.add_option('--Wjets',      dest='Wjets',      action='store_true', default=False)
    parser.add_option('--multijet',   dest='multijet',  action='store_true', default=False)
    parser.add_option('--iso',        dest='iso',       action='store_true', default=False)
    parser.add_option('--mt',         dest='mt',        action='store_true', default=False)
    parser.add_option('--data', dest='data', action='store_true', default=False)
    parser.add_option('--mc', dest='mc', action='store_true', default=False)
    return parser.parse_args()

#______________________________________________________________________________
def main():

    ops, args = options()
    assert ops.Wjets or ops.multijet or ops.iso or ops.mt, 'Specify sample.'

    ## config
    timestamp = time.strftime('%Y-%m-%d')

    ## Wenu+jets
    if ops.Wjets and ops.channel == 1:
        outfile = 'weights.tauID.el-tau.Wjets_inclusive.%s.root' % timestamp
        data_sample = data.Egamma.mydata
        mc_samples_to_subtract = [
                mc.Ztautau,
                mc.Zee,
                mc.top,
#                mc.Wgamma,
                mc.diboson,
                ]
        control_regions = control_regions_Wjets()

    ## Wmunu+jets
    elif ops.Wjets and ops.channel == 2:
        outfile = 'weights.tauID.mu-tau.Wjets_inclusive.%s.root' % timestamp
        data_sample = data.Muons.mydata
        mc_samples_to_subtract = [
                mc.Ztautau,
                mc.Zmumu,
                mc.top,
                mc.diboson,
                ]
        control_regions = control_regions_Wjets()
        
    ## multijet (e)
    elif ops.multijet and ops.channel == 1:
        outfile = 'weights.tauID.el-tau.multijet_inclusive.%s.root' % timestamp
        data_sample = data.Egamma.mydata
        mc_samples_to_subtract = [
                mc.Wenu,
                mc.Wtaunu,
                mc.Ztautau,
                mc.Zee,
                mc.top,
                mc.diboson,
                ]
        control_regions = control_regions_multijet()

    ## multijet (mu)
    elif ops.multijet and ops.channel == 2:
        outfile = 'weights.tauID.mu-tau.multijet_inclusive.%s.root' % timestamp
        data_sample = data.Muons.mydata
        mc_samples_to_subtract = [
                mc.Wmunu,
                mc.Wtaunu,
                mc.Ztautau,
                mc.Zmumu,
                mc.top,
                mc.diboson,
                ]
        control_regions = control_regions_multijet()

    ## multijet - iso (e)
    elif ops.iso and ops.channel == 1:
        outfile = 'weights.lepIso.el-tau.multijet_inclusive.%s.root' % timestamp
        data_sample = data.Egamma.mydata
        mc_samples_to_subtract = [
                mc.Wenu,
                mc.Wtaunu,
                mc.Ztautau,
                mc.Zee,
                mc.top,
                mc.diboson,
                ]
        control_regions = control_regions_iso()

    ## multijet - iso (mu)
    elif ops.iso and ops.channel == 2:
        outfile = 'weights.lepIso.mu-tau.multijet_inclusive.%s.root' % timestamp
        data_sample = data.Muons.mydata
        mc_samples_to_subtract = [
                mc.Wmunu,
                mc.Wtaunu,
                mc.Ztautau,
                mc.Zmumu,
                mc.top,
                mc.diboson,
                ]
        control_regions = control_regions_iso()

    ## W+jet - mT (e)
    elif ops.mt and ops.channel == 1:
        outfile = 'weights.tauID.el-tau.Wjet_mT.%s.root' % timestamp
        data_sample = data.Egamma.mydata
        mc_samples_to_subtract = [
                #mc.Wenu,
                #mc.Wtaunu,
                mc.Ztautau,
                mc.Zee,
                mc.top,
                mc.diboson,
                ]
        mc_samples = [
                mc.Wenu,
                mc.Wtaunu,
                ]
        control_regions = control_regions_mT()

    ## W+jet - mT  iso (mu)
    elif ops.mt and ops.channel == 2:
        outfile = 'weights.tauID.mu-tau.Wjet_mT.%s.root' % timestamp
        data_sample = data.Muons.mydata
        mc_samples_to_subtract = [
                #mc.Wmunu,
                #mc.Wtaunu,
                mc.Ztautau,
                mc.Zmumu,
                mc.top,
                mc.diboson,
                ]
        mc_samples = [
                mc.Wmunu,
                mc.Wtaunu,
                ]     
        control_regions = control_regions_mT()

        
    ## get samples
    sample_list = [data_sample] + mc_samples_to_subtract + mc_samples
    for sample in sample_list:
        sample.load(sample.search(ops.load, useName=True))

    ## make fake rate hists
    hist_names = []
    for cr in control_regions:
        hn = make_fake_factors(cr, data_sample, mc_samples_to_subtract, mc_samples, outfile)
        hist_names.append(hn)

    print '  written hists:'
    for hn in hist_names:
        print '   ', hn

    print '  Done!'

#_______________________________________________________________________________
def control_regions_Wjets():

    control_regions = [
                ( '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_mTW_os_tau1p_tauID/h2_tau_eta_vs_tau_pt',
                  '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_mTW_os_tau1p_tauDenom/h2_tau_eta_vs_tau_pt',
                  'h2_tau_eta_vs_tau_pt_os_tau1p' ),
                ( '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_mTW_ss_tau1p_tauID/h2_tau_eta_vs_tau_pt',
                  '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_mTW_ss_tau1p_tauDenom/h2_tau_eta_vs_tau_pt',
                  'h2_tau_eta_vs_tau_pt_ss_tau1p' ),
                ( '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_mTW_os_tau3p_tauID/h2_tau_eta_vs_tau_pt',
                  '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_mTW_os_tau3p_tauDenom/h2_tau_eta_vs_tau_pt',
                  'h2_tau_eta_vs_tau_pt_os_tau3p' ),
                ( '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_mTW_ss_tau3p_tauID/h2_tau_eta_vs_tau_pt',
                  '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_mTW_ss_tau3p_tauDenom/h2_tau_eta_vs_tau_pt',
                  'h2_tau_eta_vs_tau_pt_ss_tau3p' ),
                ]

    return control_regions

#_______________________________________________________________________________
def control_regions_multijet():
    control_regions = [
            ( '/regions/multijet_tauFF/lepPt26_lepIsoDenom_metLT30_mT30_lepd0GT02_os_tau1p_tauID/h2_tau_eta_vs_tau_pt',
              '/regions/multijet_tauFF/lepPt26_lepIsoDenom_metLT30_mT30_lepd0GT02_os_tau1p_tauDenom/h2_tau_eta_vs_tau_pt',
              'h2_tau_eta_vs_tau_pt_os_tau1p' ),
            ( '/regions/multijet_tauFF/lepPt26_lepIsoDenom_metLT30_mT30_lepd0GT02_ss_tau1p_tauID/h2_tau_eta_vs_tau_pt',
              '/regions/multijet_tauFF/lepPt26_lepIsoDenom_metLT30_mT30_lepd0GT02_ss_tau1p_tauDenom/h2_tau_eta_vs_tau_pt',
              'h2_tau_eta_vs_tau_pt_ss_tau1p' ),
            ( '/regions/multijet_tauFF/lepPt26_lepIsoDenom_metLT30_mT30_lepd0GT02_os_tau3p_tauID/h2_tau_eta_vs_tau_pt',
              '/regions/multijet_tauFF/lepPt26_lepIsoDenom_metLT30_mT30_lepd0GT02_os_tau3p_tauDenom/h2_tau_eta_vs_tau_pt',
              'h2_tau_eta_vs_tau_pt_os_tau3p' ),
            ( '/regions/multijet_tauFF/lepPt26_lepIsoDenom_metLT30_mT30_lepd0GT02_ss_tau3p_tauID/h2_tau_eta_vs_tau_pt',
              '/regions/multijet_tauFF/lepPt26_lepIsoDenom_metLT30_mT30_lepd0GT02_ss_tau3p_tauDenom/h2_tau_eta_vs_tau_pt',
              'h2_tau_eta_vs_tau_pt_ss_tau3p' ),
            ]
    return control_regions

#_______________________________________________________________________________
def control_regions_iso():
    control_regions = [
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_os_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_os_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_os' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_ss_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_ss_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_ss' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT01_os_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT01_os_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_lepd0GT01_os' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT01_ss_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT01_ss_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_lepd0GT01_ss' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT02_os_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT02_os_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_lepd0GT02_os' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT02_ss_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT02_ss_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_lepd0GT02_ss' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT03_os_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT03_os_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_lepd0GT03_os' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT03_ss_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT03_ss_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_lepd0GT03_ss' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT04_os_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT04_os_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_lepd0GT04_os' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT04_ss_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0GT04_ss_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_lepd0GT04_ss' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0LT02_os_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0LT02_os_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_lepd0LT02_os' ),
            ( '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0LT02_ss_lepIso/h2_lep_eta_vs_lep_pt',
              '/regions/multijet_lepIsoFF/lepPt26_ntau0_metLT30_mT30_lepd0LT02_ss_lepIsoDenom/h2_lep_eta_vs_lep_pt',
              'h2_lep_eta_vs_lep_pt_mT30_lepd0LT02_ss' ),
            ]
    return control_regions


#_______________________________________________________________________________
def control_regions_mT():

    control_regions = [
                ( '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_os_tau1p_tauID/h_mT',
                  '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_os_tau1p_tauDenom/h_mT',
                  'h_mT_os_tau1p' ),
                ( '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_ss_tau1p_tauID/h_mT',
                  '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_ss_tau1p_tauDenom/h_mT',
                  'h_mT_ss_tau1p' ),
                ( '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_os_tau3p_tauID/h_mT',
                  '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_os_tau3p_tauDenom/h_mT',
                  'h_mT_os_tau3p' ),
                ( '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_ss_tau3p_tauID/h_mT',
                  '/regions/Wjets_tauFF/lepPt26_lepIso_tauLepVeto_ss_tau3p_tauDenom/h_mT',
                  'h_mT_ss_tau3p' ),
                ]

    return control_regions

#_______________________________________________________________________________
def make_fake_factors(cr, data_sample, mc_samples_to_subtract, mc_samples, outfile):

    ops, args = options()
    h_pass_name, h_fail_name, h_fakerate_name = cr

    # get numer,denom
    if ops.data:
        h_pass = get_data(h_pass_name, data_sample, mc_samples_to_subtract)
        h_fail = get_data(h_fail_name, data_sample, mc_samples_to_subtract)
    if ops.mc:
        h_pass_mc = get_mc(h_pass_name, mc_samples)
        h_fail_mc = get_mc(h_fail_name, mc_samples)
    else:
        h_pass = get_data_minus_mc(h_pass_name, data_sample, mc_samples_to_subtract)
        h_fail = get_data_minus_mc(h_fail_name, data_sample, mc_samples_to_subtract)
    assert h_pass, 'Numer not retrieved. Exiting.'
    assert h_fail, 'Denom not retrieved. Exiting.'

    for h in (h_pass, h_fail, h_pass_mc, h_fail_mc):
        if not h.GetSumw2N():
            h.Sumw2()

#    h_ratio = h_pass.Clone('%s_ratio' % h_pass.GetName()
    h_ratio = h_pass.Clone('%s_data' % h_fakerate_name)
    divide_option = ''
    h_ratio.Divide(h_pass, h_fail, 1.0, 1.0, divide_option)
    h_ratio_mc = h_pass.Clone('%s_mc' % h_fakerate_name)
    h_ratio_mc.Divide(h_pass_mc, h_fail_mc, 1.0, 1.0, divide_option)

    metaroot.file.write(h_ratio, outfile)
    metaroot.file.write(h_ratio_mc, outfile)
    
    ## dump contents?
    if ops.dump:
        for i_y in xrange(1, h_ratio.GetNbinsY()+1):
            for i_x in xrange(1, h_ratio.GetNbinsX()+1):
                r = h_ratio.GetBinContent(i_x, i_y)
                err  = h_ratio.GetBinError(i_x, i_y)
                x_axis = h_ratio.GetXaxis()
                y_axis = h_ratio.GetYaxis()
                x_low  = x_axis.GetBinLowEdge(i_x)
                x_high  = x_axis.GetBinLowEdge(i_x+1)
                y_low  = y_axis.GetBinLowEdge(i_y)
                y_high  = y_axis.GetBinLowEdge(i_y+1)
                print 'Bin: pt=%.0f-%.0f GeV, eta=%.2f-%.2f Pass/fail: %.4f +/- %.4f' % (x_low, x_high, y_low, y_high, r, err)
        print '\n'

    return h_fakerate_name


#_______________________________________________________________________________
def get_data_minus_mc(hist_path, data_sample, mc_samples):

    ops, args = options()
    ## get data, make clone
    h = data_sample.get(hist_path)
    assert h, 'Did not retrieve data hist for %s. Exiting.' % hist_path

    if ops.mt:
        h_total = h.Clone('%s_minus_mc' % h.GetName())
    else:
        h_total = h.ProjectionX('%s_minus_mc' % h.GetName(), 1, 3)

    ## memory resident (not owned by any TFile)
    h_total.SetDirectory(0)

    ## get MC, subtract from data
    for s in mc_samples:
        h_mc = s.get(hist_path)
        if h_mc:
            ## Force each mc histogram to be positive definite (not just their sum).
            ## You should have previously merged together components that can have
            ## interference.
            if ops.mt:
                h_mc = set_negative_bins_zero(h_mc)
                h_total.Add(h_mc, -1.0)
            else:
                h_mc_proj = h_mc.ProjectionX('h_mc', 1, 3)
                h_mc_proj = set_negative_bins_zero(h_mc_proj)
                h_total.Add(h_mc_proj, -1.0)

    ## Force no negative bins in the result
    h_total = set_negative_bins_zero(h_total, warn=True)

    return h_total

#_______________________________________________________________________________
def get_data(hist_path, data_sample, mc_samples_to_subtract):

    ops, args = options()
    ## get data, make clone
    h = data_sample.get(hist_path)
    assert h, 'Did not retrieve data hist for %s. Exiting.' % hist_path

    if ops.mt:
        h_total = h.Clone('%s_minus_mc' % h.GetName())
    else:
        h_total = h.ProjectionX('%s_minus_mc' % h.GetName(), 1, 3)

    ## memory resident (not owned by any TFile)
    h_total.SetDirectory(0)

    ## get MC, subtract from data
    for s in mc_samples_to_subtract:
        h_mc = s.get(hist_path)
        if h_mc:
            ## Force each mc histogram to be positive definite (not just their sum).
            ## You should have previously merged together components that can have
            ## interference.
            if ops.mt:
                h_mc = set_negative_bins_zero(h_mc)
                h_total.Add(h_mc, -1.0)
            else:
                h_mc_proj = h_mc.ProjectionX('h_mc', 1, 3)
                h_mc_proj = set_negative_bins_zero(h_mc_proj)
                h_total.Add(h_mc_proj, -1.0)

    ## Force no negative bins in the result
    h_total = set_negative_bins_zero(h_total, warn=True)

    return h_total

#_______________________________________________________________________________
def get_mc(hist_path, mc_samples):

    ops, args = options()

    ## get MC
    first_hist = True
    
    for s in mc_samples:
        h_mc = s.get(hist_path)
        if h_mc:
            ## Force each mc histogram to be positive definite (not just their sum).
            ## You should have previously merged together components that can have
            ## interference.
            if ops.mt:
                if first_hist:
                    h_mc = set_negative_bins_zero(h_mc)
                    h_total = h_mc.Clone('%s_mc' % h_mc.GetName())
                    h_total.SetDirectory(0)
                    first_hist = False
                else:
                    h_mc = set_negative_bins_zero(h_mc)
                    h_total.Add(h_mc, 1.0)
            else:
                h_mc_proj = h_mc.ProjectionX('h_mc', 1, 3)
                h_mc_proj = set_negative_bins_zero(h_mc_proj)
                h_total.Add(h_mc_proj, -1.0)

    ## Force no negative bins in the result
    h_total = set_negative_bins_zero(h_total, warn=True)

    return h_total

#______________________________________________________________________________
def set_negative_bins_zero(h, warn=False):
    for bin in xrange(0, h.GetNbinsX()+2):
        if h.GetBinContent(bin) < 0:
            h.SetBinContent(bin, 0.0)
            h.SetBinError(bin, 0.0)  ## NOTE
            if warn:
                print 'WARNING (calc_fakefactors.py): setting negative bin to zero in'
                print '  %s' % h.GetName()
    return h


#_____________________________________________________________________________
if __name__ == '__main__': main()

