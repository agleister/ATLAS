#!/usr/bin/env python
"""
plot_with_sys.py

author: Ryan Reece <ryan.reece@cern.ch>
created: 2011-06-01
mofified from plot.py to include systematics by Andrew Leister <aleister@cern.ch>

run with plot_with_sys.py --blind --ratio --kfactors -c ${CHANNEL} --load=$CONDOR_TOP --faketaus -p "/regions/Zprime/.*" --comment="Zprime"
   where $CONDOR_TOP is the directory containing $CONDOR_MAIN and the systematic variations

"""
#------------------------------------------------------------------------------

## std
import math
import optparse
import os
import re
import time
import glob

## ROOT
import ROOT, rootlogon
ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = 1001

# turn off warnings
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

## my modules
import metaroot
from redict import redict
import fileutils
import poissonize

## samples
from samples.p1344 import mc
from samples.p1443 import data

timestamp = time.strftime('%Y-%m-%d-%Hh%Mm%Ss')

# config blinding
blind_patterns = []
varbins_blind = 395.5
blinding_dict = {
        'h_mTtot_varbins'    : varbins_blind,
        'h_mTtot_high'       : 400.,
        'h_mTtot_mid_high'   : 400.,
        'h_mTtot_integrated' : 400.,
        'h_mVis_high'        : 400.,
        'h_mVis_varbins'     : varbins_blind,
        'h_mT_high'          : 350.,
        'h_tau_pt_high'      : 250.,
        'h_tau_pt_varbins'   : 227.57,
        'h_lep_pt_high'      : 200.,
        'h_lep_pt_varbins'   : 189.3,
        'h_met_high'         : 200.,
        'h_met_varbins'      : 168.0,
        'h_mmc_mass_high'    : 700.,
        'h_mmc_mass_varbins' : varbins_blind,
        'h_col_mass_high'    : 700.,
        'h_col_mass_varbins' : varbins_blind,
        }    

#_____________________________________________________________________________
def options():
    parser = optparse.OptionParser(description="options")
    # on/off switches
    parser.add_option('-b', '--blind'  ,   dest='blind'  ,   action='store_true', default=False)
    parser.add_option('-D', '--no-data'  , dest='no_data'  , action='store_true', default=False)
    parser.add_option('-k', '--kfactors' , dest='kfactors' , action='store_true', default=False)
    parser.add_option('-r', '--ratio'    , dest='ratio'    , action='store_true', default=False)
    parser.add_option('-S', '--no-signal', dest='no_signal', action='store_true', default=False)
    # configurables
    parser.add_option('-B', '--blind-pattern-file', dest='blind_pattern_file', type=str, default=None) # does nothing right now
    parser.add_option('-c', '--channel'     , dest='channel'     , type=int, default=None)
    parser.add_option('-C', '--comment'     , dest='comment'     , type=str, default=None)
    parser.add_option('-o', '--output-file' , dest='output_file' , type=str, default=None)
    parser.add_option('-l', '--load'        , dest='load'        , type=str, default=None)
    parser.add_option('-p', '--pattern'     , dest='pattern'     , type=str, default=None)
    parser.add_option('-P', '--pattern-file', dest='pattern_file', type=str, default=None)
    parser.add_option('-s', '--signalxsec'  , dest='signalxsec'  , type=int, default=1)
    # bkg estimation
    parser.add_option('--OSminusSS', dest='OSminusSS', action='store_true', default=False)
    parser.add_option('--faketaus',  dest='faketaus',  action='store_true', default=False)
    parser.add_option('--faketaus2', dest='faketaus2', action='store_true', default=False)
    parser.add_option('--nominal',   dest='nominal',   action='store_true', default=False)
    parser.add_option('--multijet',  dest='multijet',  action='store_true', default=False)
    parser.add_option('--higgs',     dest='higgs',     action='store_true', default=False)
    return parser.parse_args()

#______________________________________________________________________________
def main():
    ops, args = options()

    ## config
    if not ops.channel in [1, 2]                          : sys.exit(' ERROR: Please provide --channel 1 (el-had) or 2 (mu-had).')
#    if not (ops.OSminusSS or ops.faketaus or ops.faketaus2 or ops.nominal) : sys.exit(' ERROR: Please specify background estimation --OSminusSS or --faketaus or --nominal.')
    if not ops.load                                       : sys.exit(' ERROR: Please provide --load path for samples.')

    if ops.channel == 1:
        channel_name = 'el-tau'
        lepton_name  = 'e'
    if ops.channel == 2:
        channel_name = 'mu-tau'
        lepton_name  = '#mu'

    project_name = 'Zp12'
    if ops.output_file:
        output_file = ops.output_file
    else:
        if ops.comment:
            output_file = '%s.%s.%s.%s.canv.root' %  (project_name, channel_name, timestamp, ops.comment)
        else:
            output_file = '%s.%s.%s.canv.root'    %  (project_name, channel_name, timestamp)

    data_sample = None
    bkg_samples = []
    #data_sample_bch = {}
    #bkg_samples_bch = {}
    #bkg_samples_bch['bch'] = []

    ## electron samples
    if ops.channel == 1:

        # data 
        data_sample = data.Egamma.mydata
        #data_sample_bch['bch'] = data.Egamma.mydata

        # backgrounds
        if ops.OSminusSS:
            bkg_samples.append( mc.OSminusSS['Ztautau'] )
            bkg_samples.append( mc.OSminusSS['Zee'    ] )
            bkg_samples.append( mc.OSminusSS['Wenu'   ] )
            bkg_samples.append( mc.OSminusSS['Wtaunu' ] )
            bkg_samples.append( mc.OSminusSS['top'    ] )
            bkg_samples.append( mc.OSminusSS['diboson'] )
            bkg_samples.append( data.SameSign )
        else:
            bkg_samples.append( mc.Ztautau )
            if ops.nominal:
                bkg_samples.append( mc.Wenu   )
                bkg_samples.append( mc.Wtaunu )
            if ops.faketaus:
                bkg_samples.append( data.WjetsFakeTaus )
            elif ops.faketaus2:
                bkg_samples.append( data.WjetsFakeTaus )
                bkg_samples.append( data.MultijetFakeTaus )
            if ops.multijet:
                bkg_samples.append( data.MultijetFakeTaus )
            bkg_samples.append( mc.Zee     )
            bkg_samples.append( mc.top     )
            bkg_samples.append( mc.diboson )
        #bkg_samples_bch = list(bkg_samples)
        #bkg_samples_bch['bch'].append( mc.Ztautau )
        #bkg_samples_bch['bch'].append( data.WjetsFakeTaus )
        #bkg_samples_bch['bch'].append( mc.Zee     )
        #bkg_samples_bch['bch'].append( mc.top     )
        #bkg_samples_bch['bch'].append( mc.diboson )

        #systematic backgrounds
        systematics_updn = ['ELE_SCALE', 'ELE_TRIG', 'JES', 'MET_ScaleSoftTerms',
                            'TAUID', 'TESREAL', 'TAU3P', 'ELE_EFF',
                            'TESFAKE', 'TAU_ISELECTRON', 'KFACTOR']
        systematics_up = ['ELE_RES_UP', 'MET_ResoSoftTerms_UP']


        bkg_samples_sys = {}
        for sys in systematics_up:
            bkg_samples_sys[sys] = []
            bkg_samples_sys[sys].append(mc.Ztautau)
            bkg_samples_sys[sys].append(mc.Zee)
            bkg_samples_sys[sys].append(mc.top)
            bkg_samples_sys[sys].append(mc.diboson)
        for sys in systematics_updn:
            for updn in ['UP', 'DN']:
                bkg_samples_sys['%s_%s' % (sys,updn)] = []
                bkg_samples_sys['%s_%s' % (sys,updn)].append(mc.Ztautau)
                bkg_samples_sys['%s_%s' % (sys,updn)].append(mc.Zee)
                bkg_samples_sys['%s_%s' % (sys,updn)].append(mc.top)
                bkg_samples_sys['%s_%s' % (sys,updn)].append(mc.diboson)
        

        # signal
        #   signal_samples_bch = {}
        if ops.no_signal:
            signal_samples = []
        else:
            if ops.higgs:
                signal_samples = [
                                    mc.ggH125,
                                    mc.VBFH125,
                                 ]
            else: # Zprime
                signal_samples = [
#                                    mc.ZprimePythia1000,
#                                    mc.ZprimePythia1250,
                                    mc.ZprimeSSM1000,
                                    #mc.ZprimeSSM1250,
                                    mc.ZprimeSSM1500,
                                    mc.ZprimeSSM2000,
                                 ]
        #   signal_samples_bch['bch'] = list(signal_samples)

    ## muon samples
    if ops.channel == 2:

        # data 
        data_sample = data.Muons.mydata
        #data_sample_bch['bch'] = data.Muons.mydata

        # backgrounds
        if ops.OSminusSS:
            bkg_samples.append( mc.OSminusSS['Ztautau'] )
            bkg_samples.append( mc.OSminusSS['Zmumu'  ] )
            bkg_samples.append( mc.OSminusSS['Wmunu'  ] )
            bkg_samples.append( mc.OSminusSS['Wtaunu' ] )
            bkg_samples.append( mc.OSminusSS['top'    ] )
            bkg_samples.append( mc.OSminusSS['diboson'] )
            bkg_samples.append( data.SameSign )
        else:
            bkg_samples.append( mc.Ztautau )
            if ops.nominal:
                bkg_samples.append( mc.Wmunu  )
                bkg_samples.append( mc.Wtaunu )
            if ops.faketaus:
                bkg_samples.append( data.WjetsFakeTaus )
            elif ops.faketaus2:
                bkg_samples.append( data.WjetsFakeTaus )
                bkg_samples.append( data.MultijetFakeTaus )
            if ops.multijet:
                bkg_samples.append( data.MultijetFakeTaus )
            bkg_samples.append( mc.Zmumu   )
            bkg_samples.append( mc.top     )
            bkg_samples.append( mc.diboson )
        #bkg_samples_bch = list(bkg_samples)
        #bkg_samples_bch['bch'].append( mc.Ztautau )
        #bkg_samples_bch['bch'].append( data.WjetsFakeTaus )
        #bkg_samples_bch['bch'].append( mc.Zmumu     )
        #bkg_samples_bch['bch'].append( mc.top     )
        #bkg_samples_bch['bch'].append( mc.diboson )

        #systematic backgrounds
        systematics_updn = ['MUON_SCALE', 'TESREAL', 'TESFAKE', 'JES',
                            'MET_ScaleSoftTerms', 'MU_EFF', 'TAU3P', 'MU_TRIG',
                            'TAU_ISMUON', 'KFACTOR', 'TAUID',]
        systematics_up = ['MUON_ID_UP', 'MUON_MS_UP', 'MET_ResoSoftTerms_UP'] 

        bkg_samples_sys = {}
        for sys in systematics_up:
            bkg_samples_sys[sys] = []
            bkg_samples_sys[sys].append(mc.Ztautau)
            bkg_samples_sys[sys].append(mc.Zmumu)
            bkg_samples_sys[sys].append(mc.top)
            bkg_samples_sys[sys].append(mc.diboson)
        for sys in systematics_updn:
            for updn in ['UP', 'DN']:
                bkg_samples_sys['%s_%s' % (sys,updn)] = []
                bkg_samples_sys['%s_%s' % (sys,updn)].append(mc.Ztautau)
                bkg_samples_sys['%s_%s' % (sys,updn)].append(mc.Zmumu)
                bkg_samples_sys['%s_%s' % (sys,updn)].append(mc.top)
                bkg_samples_sys['%s_%s' % (sys,updn)].append(mc.diboson)
            

        # signal
        #   signal_samples_bch = {}
        if ops.no_signal:
            signal_samples = []
        else:
            if ops.higgs:
                signal_samples = [
                                    mc.ggH125,
                                    mc.VBFH125,
                                 ]
            else: # Zprime
                signal_samples = [
#                                    mc.ZprimePythia1000,
#                                    mc.ZprimePythia1250,
                                    mc.ZprimeSSM1000,
                                    #mc.ZprimeSSM1250,
                                    mc.ZprimeSSM1500,
                                    mc.ZprimeSSM2000,
                                 ]
        #   signal_samples_bch['bch'] = list(signal_samples)

    ## scale signal
    if ops.signalxsec != 1:
        sys.exit(' ERROR: Not configured for signal xsec >1.')
        #for signal_sample in signal_samples:
        #    signal_sample.latex += '(%ixSM)' % ops.signalxsec
        #    for sub_sample in signal_sample.sub_samples:
        #        sub_sample.cross_section*=float(ops.signalxsec)

    ## load samples
    #print
    #print ' @@@@@@ Loading @@@@@@'
    #for sample in [data_sample] + bkg_samples + signal_samples:
    #    print '', sample.name
    #    sample.load(sample.search(ops.load, useName=True))
    #print

    ## load main samples 
    print
    print ' @@@@@@ Loading @@@@@@'
    for sample in [data_sample]: #+ bkg_samples + signal_samples:
        print '', sample.name
        print ops.load
        print glob.glob('%s/*MAIN*' % ops.load)
        print "{run}/merged/".format(run=glob.glob('%s/*MAIN*' % ops.load)[0])
        sample.load(sample.search('{0}/merged/'.format(glob.glob('%s/*MAIN*' % ops.load)[0]), useName=True))
        #print glob.glob('%s/*BCH_CLEANING*' % ops.load)
        #print "{run}/merged/".format(run=glob.glob('%s/*BCH_CLEANING*' % ops.load)[0])
        #sample.load(sample.search('{0}/merged/'.format(glob.glob('%s/*BCH_CLEANING*' % ops.load)[0]), useName=True))
    print

    ## load bch-cleaned samples 
    #print
    #print ' @@@@@@ Loading @@@@@@'
    #for sample in [data_sample_bch['bch']]: #+ bkg_samples_bch['bch'] + signal_samples_bch['bch']:
    #    print '', sample.name
    #    print glob.glob('%s/*BCH_CLEANING*' % ops.load)
    #    print "{run}/merged/".format(run=glob.glob('%s/*BCH_CLEANING*' % ops.load)[0])
    #    sample.load(sample.search('{0}/merged/'.format(glob.glob('%s/*BCH_CLEANING*' % ops.load)[0]), useName=True))
    #print

    ## load systematics samples
    #print
    #print ' @@@@@@ Loading @@@@@@'
    #for sys in systematics_up:
    #    for sample in bkg_samples_sys[sys]:
    #        print '', sample.name
    #        sample.load(sample.search('{0}/merged'.format(glob.glob('%s/*%s*' % (ops.load, sys))[0]), useName=True))
    #for sys in systematics_updn:
    #    for updn in ['UP', 'DN']:
    #        for sample in bkg_samples_sys['%s_%s' % (sys, updn)]:
    #            print '', sample.name
    #            sample.load(sample.search('{0}/merged'.format(glob.glob('%s/*%s_%s*' % (ops.load, sys, updn))[0]), useName=True))
    #print


    ## get plot patterns
    patterns = []
    if ops.pattern:
        print ' @@@@@@ Found pattern: %s @@@@@@' % ops.pattern
        print
        patterns.append(ops.pattern)
    if ops.pattern_file:
        f = open(ops.pattern_file)
        for line in f:
            line = line.split('#')[0].strip()
            if not line: continue
            patterns.append(line)
        f.close()
        print ' @@@@@@ Found pattern file (%s) with %s patterns. @@@@@@' % (ops.pattern_file, len(patterns))
        print patterns
        print
        
    ## get blind patterns
    if ops.blind_pattern_file:
        f = open(ops.blind_pattern_file)
        for line in f:
            line = line.split('#')[0].strip()
            if not line: continue
            blind_patterns.append(line)
        f.close()
        print ' @@@@@@ Found blind pattern file (%s) with %s patterns. @@@@@@' % (ops.blind_pattern_file, len(blind_patterns))
        print
        
    t_start = time.time()
    n_plots = 0

    ## walk directories and plot canvases
    print ' @@@@@@ Walking @@@@@@'
    for dirpath, dirnames, objnames in data_sample.walk():
        for name in objnames:
            hist_path = os.path.join(dirpath, name)
            if patterns:
                if not any(re.match(p, hist_path) for p in patterns):
                    continue
            print hist_path
            plot_method = plotting_methods.rule(hist_path)
            if isinstance(plot_method, list):
                plot_methods_todo = plot_method
            else:
                 plot_methods_todo = [plot_method]
            for pm in plot_methods_todo:
                #if ("event" in hist_path and "mTtot" in hist_path and not "count" in hist_path) or () or () or () :
                #if not "count" in hist_path:
                if ("event" in hist_path and not "count" in hist_path) or ("/lep/" in hist_path) or ("/tau/" in hist_path) or ("/met" in hist_path
                     and not "along" in hist_path) :
                #if "numTrack" in hist_path:
                    if pm:
                        pm(ops, hist_path, data_sample, bkg_samples, signal_samples, #data_sample_bch, bkg_samples_bch, signal_samples_bch,
                           bkg_samples_sys, systematics_up, systematics_updn, output_file)
                        n_plots += 1
            if n_plots % 100 == 0:
                print ' %4s plots' % n_plots
                
    ## close samples
    sample_list = [data_sample] + bkg_samples + signal_samples
    #sample_list += [data_sample_bch] + bkg_samples_bch + signal_samples_bch
    #for sys in systematics_up:
    #    sample_list += bkg_samples_sys[sys]
    #for sys in systematics_updn:
    #    for updn in ['UP', 'DN']:
    #        sample_list += bkg_samples_sys['%s_%s' % (sys, updn)]
    for sample in sample_list:
        print "closing %s" % sample
        sample.close()
    t_stop = time.time()
    print
    print ' @@@@@@ Summary @@@@@@'
    print '  # plots    = %i' % n_plots
    print '  time spent = %i s' % (t_stop-t_start)
    print '  avg rate   = %.2f Hz' % (float(n_plots)/(t_stop-t_start))
    print
    print ' Done! ^.^'
    print

            
#-------------------------------------------------------------------------------
# plotting methods
#-------------------------------------------------------------------------------
ops, args = options()

plotting_methods = redict()

#_______________________________________________________________________________
def plot_simple(ops, hist_path,
                data_sample, bkg_samples, signal_samples, 
                #data_sample_bch, bkg_samples_bch, signal_samples_bch, 
                bkg_samples_sys, 
                systematics_up, systematics_updn,
                output_file,
                log_y = False,
                log_x = False
                ):
    ops, args = options()

    hist_dir = os.path.dirname(hist_path)
    #print "hist_dir = %s" % hist_dir
    hist_name = os.path.basename(hist_path)
    #print "hist_name = %s" % hist_name
    ## hack for misplacement of dphi hists
    hist_path_new = hist_path
    #if hist_name == "h_lepmet_dphi" or hist_name == "h_taumet_dphi":
    #    hist_path = "%sevent/%s" % (hist_dir[:-3], hist_name)
    #    print "original hist path = %s" % hist_path_new
    #    print "adjusted hist_path = %s" % hist_path
    data_hist = data_sample.get(hist_path_new, rename=data_sample.name)
    print data_sample.name

    #load main background hists
    print
    print ' @@@@@@ Loading main @@@@@@'
    for sample in bkg_samples + signal_samples:
        print '', sample.name
        print ops.load
        #print glob.glob('%s/*MAIN*' % ops.load)
        #print "{run}/merged/".format(run=glob.glob('%s/*MAIN*' % ops.load)[0])
        sample.load(sample.search('{0}/merged/'.format(glob.glob('%s/*MAIN*' % ops.load)[0]), useName=True))
    print
    bkg_hists,    bkg_plot_options,    bkg_labels    = get_hists_and_options(hist_path, bkg_samples   )
    print "%i bg samples" % len(bkg_samples)
    print "%i bg hists" % len(bkg_hists)
    signal_hists, signal_plot_options, signal_labels = get_hists_and_options(hist_path, signal_samples)
    #for sample in bkg_samples:
    #    sample.close()
    
    #load bch-cleaned background and signal hists
    #print
    #print ' @@@@@@ Loading bch @@@@@@'
    #for sample in bkg_samples_bch['bch'] + signal_samples_bch['bch']:
    #    print '', sample.name
        #print glob.glob('%s/*BCH_CLEANING*' % ops.load)
        #print "{run}/merged/".format(run=glob.glob('%s/*BCH_CLEANING*' % ops.load)[0])
    #    sample.load(sample.search('{0}/merged/'.format(glob.glob('%s/*BCH_CLEANING*' % ops.load)[0]), useName=True))
    #print
    #data_hist_bch = data_sample_bch['bch'].get(hist_path_new, rename=data_sample_bch['bch'].name)
    #print hist_name
    #data_hist_bch = data_hist.Clone("data_hist_bch")
    #bkg_hists_bch,    bkg_plot_options_bch,    bkg_labels_bch    = get_hists_and_options(hist_path_new, bkg_samples_bch['bch']   )
    #signal_hists_bch, signal_plot_options_bch, signal_labels_bch = get_hists_and_options(hist_path_new, signal_samples_bch['bch'])
    #for sample in bkg_samples_bch['bch'] + signal_samples_bch['bch']:
    #    sample.close()

    #load systematics hists
    bkg_hists_sys = {}
    bkg_plot_options_sys = {}
    bkg_labels_sys = {}
    for sys in systematics_up:
        for sample in bkg_samples_sys[sys]:
            sample.load(sample.search('{0}/merged'.format(glob.glob('%s/*%s*' % (ops.load, sys))[0]), useName=True))
        bkg_hists_sys[sys], bkg_plot_options_sys[sys], bkg_labels_sys[sys] = get_hists_and_options(hist_path, bkg_samples_sys[sys])
    for sys in systematics_updn:
        for updn in ['UP', 'DN']:
            for sample in bkg_samples_sys['%s_%s' % (sys, updn)]:
                sample.load(sample.search('{0}/merged'.format(glob.glob('%s/*%s_%s*' % (ops.load, sys, updn))[0]), useName=True))
            bkg_hists_sys['%s_%s' % (sys, updn)], bkg_plot_options_sys['%s_%s' % (sys, updn)], bkg_labels_sys['%s_%s' % (sys, updn)] = get_hists_and_options(hist_path, bkg_samples_sys['%s_%s' % (sys, updn)])
            

    canvas_options = metaroot.hist.CanvasOptions(log_x=log_x, log_y=log_y)

    ## config blinding
    no_data = ops.no_data
    if not no_data and data_hist and ops.blind:
        blind_data_hist = blind_the_data(data_hist, hist_name)
        data_hist = blind_data_hist

     ## config blinding for bch 
    #no_data = ops.no_data
    #if not no_data and data_hist_bch and ops.blind:
    #    blind_data_hist = blind_the_data(data_hist_bch, hist_name)
    #    data_hist_bch = blind_data_hist

#    if not no_data and data_hist and blind_patterns \
#            and any(re.match(p, hist_path) for p in blind_patterns):
#        blind_data_hist = blind_the_data(data_hist, hist_name)
#        data_hist = blind_data_hist

    ## rebinning
    coarseness = 2.0 # if not 'VBF' in hist_path else 5.0
    data_hist    = do_the_great_rebinning([data_hist] , hist_name, coarseness, log_x=log_x, log_y=log_y)[0]
    bkg_hists    = do_the_great_rebinning(bkg_hists   , hist_name, coarseness, log_x=log_x, log_y=log_y)
    signal_hists = do_the_great_rebinning(signal_hists, hist_name, coarseness, log_x=log_x, log_y=log_y)
    all_hists = [data_hist] + bkg_hists + signal_hists
    #all_hists = [data_hist] + bkg_hists #+ signal_hists
    #data_hist_bch    = do_the_great_rebinning([data_hist_bch] , hist_name, coarseness, log_x=log_x, log_y=log_y)[0]
    #bkg_hists_bch    = do_the_great_rebinning(bkg_hists_bch   , hist_name, coarseness, log_x=log_x, log_y=log_y)
    #signal_hists_bch = do_the_great_rebinning(signal_hists_bch, hist_name, coarseness, log_x=log_x, log_y=log_y)
    #all_hists += [data_hist_bch] + bkg_hists_bch + signal_hists_bch
    #all_hists +=  bkg_hists_bch + signal_hists_bch
                    
    for sys in systematics_up:
        bkg_hists_sys[sys] = do_the_great_rebinning(bkg_hists_sys[sys]   , hist_name, coarseness, log_x=log_x, log_y=log_y)
        all_hists += bkg_hists_sys[sys]
    for sys in systematics_updn:
        for updn in ['UP', 'DN']:
            bkg_hists_sys['%s_%s' % (sys, updn)] = do_the_great_rebinning(bkg_hists_sys['%s_%s' % (sys, updn)]   , hist_name, coarseness, log_x=log_x, log_y=log_y)
            all_hists += bkg_hists_sys['%s_%s' % (sys, updn)]

    ## remove negative valued bins
    for h in all_hists:
        for i_x in xrange(0, h.GetNbinsX()+2):
            if h.GetBinContent(i_x) < 0:
#                print 'WARNING (plot.py): setting negative bin to zero in'
#                print '  %s' % hist_path
#                print '  %s' % h.GetName()
                h.SetBinContent(i_x, 0.0)
                h.SetBinError(i_x, 0.0)

    ## ranges, titles
    y_max = metaroot.hist.set_max(all_hists)
    if log_y:
        y_max *= 10.0
    #title_x = data_hist.GetXaxis().GetTitle()
    #title_y = data_hist.GetYaxis().GetTitle()
    title_x = data_hist.GetXaxis().GetTitle()
    title_y = data_hist.GetYaxis().GetTitle()

    ## plot mc stack

    ###before bch cleaning
    plot = metaroot.hist.stack_hists(
            hists = bkg_hists,
            name = hist_name,
            title = ';%s;%s' % (title_x, title_y),
            min = 8e-2 if log_y else 0.0,
            max = y_max,
            plot_options = bkg_plot_options,
            canvas_options = canvas_options,
            )
    
    ###after bch cleaning
    #plot = metaroot.hist.stack_hists(
    #        hists = bkg_hists_bch,
    #        name = hist_name,
    #        title = ';%s;%s' % (title_x, title_y),
    #        min = 8e-2 if log_y else 0.0,
    #        max = y_max,
    #        plot_options = bkg_plot_options,
    #        canvas_options = canvas_options,
    #        )

    ## calc stat error
    h_error = None
    if bkg_hists:
        h_stat_error = bkg_hists[0].Clone("h_stat_error")
        for bh in bkg_hists[1:]:
            h_stat_error.Add(bh)

        ## build error hists
        error_h_list = [h_stat_error]
        #setattr(h_stat_error, 'my_plot_options', metaroot.hist.PlotOptions(fill_color=ROOT.kOrange, line_width=0, marker_size=0.0))
        h_error = h_stat_error.Clone('h_error')
        set_hist_error_in_quadrature(h_error, error_h_list)

        ## draw error
        #h_error.SetFillStyle(3254)
        #h_error.SetFillColor(ROOT.kGray+3)
        #h_error.SetMarkerStyle(0)
        #h_error.Draw("SAME,E2")

    elif data_hist:
        h_error = data_hist.Clone('h_error')
        h_error.Reset()

    ## calc scaling factors for bch error bands
    #h_bch_stack = None
    #if bkg_hists_bch:
    #    h_bch_stack = bkg_hists_bch[0].Clone("h_bch_error")
    #    for bh in bkg_hists_bch[1:]:
    #        h_bch_stack.Add(bh)
    #    h_bch_scale = h_bch_stack.Clone("h_bch_scale")
    #    for ibin in xrange(1, h_stat_error.GetNbinsX()+1):
    #        if not h_stat_error.GetBinContent(ibin) > 0.0:
    #            h_bch_scale.SetBinContent(ibin, 0.0)
    #        else:
    #            h_bch_scale.SetBinContent(ibin, h_bch_stack.GetBinContent(ibin) / h_stat_error.GetBinContent(ibin))

    ## calc sys error
    if bkg_hists:
        h_sys_error = h_error.Clone("h_sys_error")
        h_tot_error = h_error.Clone("h_tot_error")
        #   h_tot_error = h_bch_stack.Clone("h_tot_error")
        #Get sum for systematic varied hists:
        bkg_histsum_sys = {}
        fake_index = bkg_samples.index(data.WjetsFakeTaus)
        print 'fake index = %i' % fake_index 
        Wjet_plot = bkg_hists[fake_index].Clone("h_Wjet_sys")
        for sys in systematics_updn:
            bkg_histsum_sys['%s_UP' % sys] = bkg_hists_sys['%s_UP' % sys][0].Clone("h_sys_error_%s_UP" % sys)
            bkg_histsum_sys['%s_DN' % sys] = bkg_hists_sys['%s_DN' % sys][0].Clone("h_sys_error_%s_DN" % sys)
            for bh in bkg_hists_sys['%s_UP' % sys][1:]:
                bkg_histsum_sys['%s_UP' % sys].Add(bh)
            bkg_histsum_sys['%s_UP' % sys].Add(Wjet_plot)
            for bh in bkg_hists_sys['%s_DN' % sys][1:]:
                bkg_histsum_sys['%s_DN' % sys].Add(bh)
            bkg_histsum_sys['%s_DN' % sys].Add(Wjet_plot)
        for sys in systematics_up:
            bkg_histsum_sys[sys] = bkg_hists_sys[sys][0].Clone("h_sys_error_%s" % sys)
            for bh in bkg_hists_sys[sys][1:]:
                bkg_histsum_sys[sys].Add(bh)
            bkg_histsum_sys[sys].Add(Wjet_plot)
        
        #calculate up and down errors   
        for ibin in xrange(1, h_sys_error.GetNbinsX()+1):
            n = h_sys_error.GetBinContent(ibin)
            #en = h_sys_error.GetBinError(ibin)
            #esys_sum2 = en**2
            esys_sum2 = 0.0
            for sys in systematics_updn:
                n_up = bkg_histsum_sys['%s_UP' % sys].GetBinContent(ibin)
                n_dn = bkg_histsum_sys['%s_DN' % sys].GetBinContent(ibin)
                esys = (abs(n_up - n) + abs(n_dn - n)) / 2.0
                esys_sum2 += esys**2
            for sys in systematics_up:
                n_up = bkg_histsum_sys[sys].GetBinContent(ibin)
                esys = abs(n_up - n)
                esys_sum2 += esys**2
            #add the 36% uncertainty from W+jets
            esys = 0.36 * Wjet_plot.GetBinContent(ibin)
            esys_sum2 += esys**2
            h_sys_error.SetBinError(ibin, math.sqrt(esys_sum2))
            h_tot_error.SetBinError(ibin, math.sqrt(esys_sum2 + h_error.GetBinError(ibin)**2))
            #h_tot_error.SetBinError(ibin, h_bch_scale.GetBinContent(ibin)*math.sqrt(esys_sum2 + h_error.GetBinError(ibin)**2))

        ## build error hists
        error_h_list = []
        error_h_list.extend( [h_stat_error, h_sys_error])
        #error_h_list.extend( [h_sys_error])
        setattr(h_sys_error, 'my_plot_options', metaroot.hist.PlotOptions(fill_color=ROOT.kYellow, line_width=0, marker_size=0.0))
        setattr(h_stat_error, 'my_plot_options', metaroot.hist.PlotOptions(fill_color=ROOT.kOrange+10, line_width=0, marker_size=0.0))
        h_error = h_tot_error.Clone('h_error')
        #h_nobch_error = h_sys_error.Clone("h_nobch_error")
        #set_hist_error_in_quadrature(h_error, error_h_list)

        ## draw error
        h_error.SetFillStyle(3254)
        h_error.SetFillColor(ROOT.kGray+3)
        h_error.SetMarkerStyle(0)
        h_error.Draw("SAME,E2")

    ## decide coordinates
    if is_left_sided(data_hist):
    #if is_left_sided(data_hist_bch):
        ## right side legend
        legend_x2 = 0.92
        legend_y2 = 0.88
        signal_legend_x2 = legend_x2 - 0.18
        signal_legend_y2 = legend_y2 - 0.05
        lumi_text_x = signal_legend_x2 - 0.20
        lumi_text_y = signal_legend_y2 + 0.03
        atlas_text_x = 0.20
        atlas_text_y = lumi_text_y
    else:
        ## left side legend
        legend_x2 = 0.36
        legend_y2 = 0.88
        signal_legend_x2 = legend_x2 + 0.18
        signal_legend_y2 = legend_y2 - 0.05
        lumi_text_x = signal_legend_x2 - 0.16
        lumi_text_y = signal_legend_y2 + 0.03
        atlas_text_x = 0.65
        atlas_text_y = lumi_text_y

    ## draw legend
    if not data_hist or no_data:
        legend_hists = bkg_hists + [h_error]
#        legend_labels = bkg_labels + ['stat. #oplus syst.']
        legend_labels = bkg_labels + ['uncertainty']
        legend_draw_options = ['HIST']*len(bkg_hists) + ['HIST']
    else:
        legend_hists = [data_hist] + bkg_hists + [h_error]
#        legend_labels = ['data 2012'] + bkg_labels + ['stat. #oplus syst.']
        legend_labels = ['data 2012'] + bkg_labels + ['uncertainty']
        legend_draw_options = ['PE'] + ['HIST']*len(bkg_hists) + ['HIST']
    legend = metaroot.hist.make_legend(
             hists = legend_hists,
             labels = legend_labels,
             draw_options = legend_draw_options,
             width=0.17, 
             height=0.05,
             x2=legend_x2,
             y2=legend_y2)
    legend.Draw()

    ## draw legend (bch version)
    #if not data_hist or no_data:
    #    legend_hists = bkg_hists_bch + [h_error]
#        legend_labels = bkg_labels + ['stat. #oplus syst.']
    #    legend_labels = bkg_labels + ['stat.']
    #    legend_draw_options = ['HIST']*len(bkg_hists) + ['HIST']
    #else:
    #    legend_hists = [data_hist] + bkg_hists_bch + [h_error]
#        legend_labels = ['data 2012'] + bkg_labels + ['stat. #oplus syst.']
    #    legend_labels = ['data 2012'] + bkg_labels_bch + ['uncertainty']
    #    legend_draw_options = ['PE'] + ['HIST']*len(bkg_hists_bch) + ['HIST']
        #print "%i legend_hists" % len(legend_hists)
        #print "    %i bkg_hists_bch" % len(bkg_hists_bch)
        #print "%i legend_labels" % len(legend_labels)
        #print "    %i bkg_labels" % len(bkg_labels)
        #for label in bkg_labels:
        #    print label
    #legend = metaroot.hist.make_legend(
    #        hists = legend_hists,
    #        labels = legend_labels,
    #        draw_options = legend_draw_options,
    #        width=0.17, 
    #        height=0.05,
    #        x2=legend_x2,
    #        y2=legend_y2)
    #legend.Draw()

    ## draw signal and legend
    if not ops.no_signal:
        for h, po in zip(signal_hists, signal_plot_options):
            po.configure(h)
            h.Draw('HIST same')

        signal_legend = metaroot.hist.make_legend(
                hists = signal_hists,
                labels = signal_labels,
                draw_options = ['HIST']*len(signal_hists),
                width=0.15, 
                height=0.05,
                x2=signal_legend_x2,
                y2=signal_legend_y2)
        signal_legend.Draw()

    ## draw signal and legend (bch used)
    #if not ops.no_signal:
    #    for h, po in zip(signal_hists_bch, signal_plot_options_bch):
    #        po.configure(h)
    #        h.Draw('HIST same')

    #    signal_legend = metaroot.hist.make_legend(
    #            hists = signal_hists_bch,
    #            labels = signal_labels_bch,
    #            draw_options = ['HIST']*len(signal_hists_bch),
    #            width=0.15, 
    #            height=0.05,
    #            x2=signal_legend_x2,
    #            y2=signal_legend_y2)
    #    signal_legend.Draw()

    ## draw watermark
    if True:
        atlas_text = metaroot.utils.make_text(x=atlas_text_x, y=atlas_text_y, text='ATLAS Internal', font=72)
        atlas_text.Draw()

    ## draw luminosity
    lumi = data_sample.int_lumi
    #lumi_text = metaroot.utils.make_lumi_text(x=lumi_text_x, y=lumi_text_y, lumi='%.1f fb^{-1}' % (lumi/1000.0), size=0.04)
    lumi_text = metaroot.utils.make_lumi_text(x=lumi_text_x, y=lumi_text_y, lumi='%.1f fb^{-1}' % (20300.0/1000.0), size=0.04)
    lumi_text.Draw()

    ## plot data
    if data_hist and not no_data:
        data_plot_options = metaroot.hist.PlotOptions()
        data_plot_options.configure(data_hist)
        data_hist.GetXaxis().SetTitle('')
        data_hist.GetYaxis().SetTitle('')
        data_hist.Draw('PE same')
        #g_data = poissonize.GetPoissonizedGraph(data_hist)
        #g_data.Draw("PE same")

    ## plot data (bch applied)
    #if data_hist and not no_data:
    #    data_plot_options = metaroot.hist.PlotOptions()
    #    data_plot_options.configure(data_hist)
    #    data_hist.GetXaxis().SetTitle('')
    #    data_hist.GetYaxis().SetTitle('')
    #    data_hist.Draw('PE same')

    plot['canvas'].Update()

    ## draw blind line
    if ops.blind and is_behind_blind(hist_name):
        blind_line = draw_blind_line(hist_name, plot['canvas'])

    if ops.ratio:
        ## top canvas
        top_plot = plot
        canv_top = top_plot['canvas']
        ## make data/(SM) ratio
        #g_ratio = poissonize.GetPoissonizedRatioGraph(data_hist, top_plot['stack'].GetStack().Last())
        g_ratio = poissonize.GetPoissonizedRatioGraph(data_hist, top_plot['stack'].GetStack().Last())
        ratio_systs_hists = build_ratio_systs_hists(h_error, error_h_list)
        error_plot_options = [ h.my_plot_options for h in error_h_list ]
        ratio_systs_plot = metaroot.hist.pile_hists(
                hists = [g_ratio] + ratio_systs_hists,
                name = hist_name+'_ratio_systs_canv',
                title = ';%s;obs. / exp.' % top_plot['stack'].GetXaxis().GetTitle(),
                draw_options = ['PE'] + len(ratio_systs_hists)*['E2'],
                plot_options = [metaroot.hist.PlotOptions()] + error_plot_options,
                min = 0.0,
                max = 2.0,
                include_y = [1.0],
                )
        ## bottom canvas
        canv_bottom = ratio_systs_plot['canvas']
        if log_x:
            canv_bottom.SetLogx(True)
        line08 = metaroot.utils.draw_horiz_line(canv_bottom, 0.8,color=ROOT.kGray+0, width=1, style=7)
        line10 = metaroot.utils.draw_horiz_line(canv_bottom, 1.0,color=ROOT.kGray+2, width=1, style=7)
        line12 = metaroot.utils.draw_horiz_line(canv_bottom, 1.2,color=ROOT.kGray+0, width=1, style=7)
        ## shared canvas
        shared_plot = metaroot.plot.plot_shared_axis(
            canv_top,
            canv_bottom,
            name = hist_name +'_ratio',
            split=0.35,
            axissep=0.02,
            ndivs = [505,503],
            canvas_options = metaroot.hist.CanvasOptions(width=900, height=900)
            )
        plot = shared_plot
        #plot = ratio_systs_plot

    ## append name with logx/y
    if canvas_options.log_y:
        newname = plot['canvas'].GetName()
        newname = fileutils.strip_root_ext(newname, exts=['_logy', '_logx', '_logxy'])
        if canvas_options.log_x and canvas_options.log_y:
            newname += '_logxy'
        elif canvas_options.log_x:
            newname += '_logx'
        elif canvas_options.log_y:
            newname += '_logy'
        plot['canvas'].SetName(newname)

    ## save plot
    metaroot.file.write(plot['canvas'], output_file, hist_dir) 


plotting_methods['.*/h_'] = plot_simple


#_______________________________________________________________________________
def plot_simple_log_y(ops, hist_path,
                      data_sample, bkg_samples, signal_samples, 
                      #data_sample_bch, bkg_samples_bch, signal_samples_bch,
                      bkg_samples_sys,
                      systematics_up, systematics_updn,
                      output_file):

    return plot_simple(ops, hist_path,
                       data_sample, bkg_samples, signal_samples, 
                       #data_sample_bch, bkg_samples_bch, signal_samples_bch,
                       bkg_samples_sys,
                       systematics_up, systematics_updn,
                       output_file,
                       log_y = True)

plotting_methods['.*/h_.*etcone.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*Etcone.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*ptcone.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*nucone.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*mu_d0_exPV.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*mu_rel_pt_diff.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*mu_z0_exPV.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*mu_.*chi2.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*el_trackd0pv.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_taulep.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_lepmet.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_dphisum.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_sumcosdphi.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*charge.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*high.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*dphi.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*Veto.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*TRT.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*EMScale.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*EMFraction.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_tau_olr*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_tau_lep_dz0'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_.*Veto.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_n_jets.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_n_bjets.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_lep_d0.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_lep_z0.*'] = [plot_simple, plot_simple_log_y]
plotting_methods['.*/h_mTtot_integrated.*'] = [plot_simple, plot_simple_log_y]

#_______________________________________________________________________________
def plot_simple_log_xy(ops, hist_path,
                       data_sample, bkg_samples, signal_samples, 
                       #data_sample_bch, bkg_samples_bch, signal_samples_bch,
                       bkg_samples_sys,
                       systematics_up, systematics_updn,
                       output_file
                       ):

    return plot_simple(ops, hist_path,
                       data_sample, bkg_samples, signal_samples, 
                       #data_sample_bch, bkg_samples_bch, signal_samples_bch,
                       bkg_samples_sys,
                       systematics_up, systematics_updn,
                       output_file,
                       log_y = True,
                       log_x = True)

#plotting_methods['.*/h_el_pt_high'] = [plot_simple, plot_simple_log_y, plot_simple_log_xy]
#plotting_methods['.*/h_tau_pt_high'] = [plot_simple, plot_simple_log_y, plot_simple_log_xy]
#plotting_methods['.*/h_met_high'] = [plot_simple, plot_simple_log_y, plot_simple_log_xy]
#plotting_methods['.*/h_trans_eff_mass_high'] = [plot_simple, plot_simple_log_y, plot_simple_log_xy]
plotting_methods['.*/h_.*varbins.*'] = [plot_simple_log_xy]

#_______________________________________________________________________________
def plot_h_cut_flow(ops, hist_path,
                    data_sample, bkg_samples, signal_samples, 
                    #data_sample_bch, bkg_samples_bch, signal_samples_bch,
                    bkg_samples_sys,
                    systematics_up, systematics_updn,
                    output_file):

    if hist_path.count('/h_cut_flow_raw'):
        lumi = 0 # means don't scale any histograms

    hist_dir = os.path.dirname(hist_path)
    hist_name = os.path.basename(hist_path)
    data_hist = data_sample.get(hist_path, rename=data_sample.name)
    bkg_hists, bkg_plot_options, bkg_labels = get_hists_and_options(hist_path, bkg_samples)
    signal_hists, signal_plot_options, signal_labels = get_hists_and_options(hist_path, signal_samples)
    
    ## ASCII table of cut flow
    tab = [ ['', 'data'] + bkg_labels ]
    hists = [data_hist] + bkg_hists
    for bin in xrange(1, data_hist.GetNbinsX()+1):
        cut_name = data_hist.GetXaxis().GetBinLabel(bin)
        row = [ cut_name ]
        for h in hists:
            row.append( h.GetBinContent(bin) )
        tab.append(row)

    import ascii_table
    ascii_table.write(tab, '%s.%s.table.txt' % (fileutils.strip_root_ext(output_file), hist_path.strip('/').replace('/', '_')))

    ## plot simple
    plot_simple_log_y(hist_path,
                      data_sample, bkg_samples, signal_samples, bkg_samples_sys,
                      systematics_up, systematics_updn,
                      output_file)
            
plotting_methods['/h_cut_flow'] = plot_h_cut_flow
plotting_methods['/h_cut_flow_raw'] = plot_h_cut_flow
plotting_methods['/h_n_events'] = plot_h_cut_flow
plotting_methods['.*truth.*'] = plot_h_cut_flow

### ignore
plotting_methods['.*/h2_'] = None
plotting_methods['.*/.*_check.*'] = None


#-------------------------------------------------------------------------------
# free functions
#-------------------------------------------------------------------------------

#_______________________________________________________________________________
def get_plot_option(sample):
    _color = {}
    # data
    # bkg
    _color['Zmumu']      = metaroot.style.green
    _color['Zee']        = metaroot.style.green
    _color['Ztautau']    = metaroot.style.light_blue
    _color['ZDYtautau']  = metaroot.style.light_blue
    _color['Wmunu']      = metaroot.style.red
    _color['Wenu']       = metaroot.style.red
    _color['Wtaunu']     = metaroot.style.orange
    _color['top']        = metaroot.style.brown
    _color['diboson']    = metaroot.style.beige
    _color['FakeTaus']   = metaroot.style.red
    _color['WjetsFakeTaus']   = metaroot.style.red
    _color['MultijetFakeTaus']   = metaroot.style.light_gray
    _color['SameSign']   = metaroot.style.yellow
    # sig
    _color['Zprime500tautau']  = metaroot.style.magenta
    _color['Zprime750tautau']  = metaroot.style.bright_blue
    _color['Zprime1000tautau'] = metaroot.style.bright_green
    _color['Zprime1250tautau'] = metaroot.style.bright_red
    _color['ZprimeSSM1000'] = metaroot.style.bright_green
    _color['ZprimeSSM1250'] = metaroot.style.bright_red
    _color['ZprimeSSM1500'] = metaroot.style.bright_red
    _color['ZprimeSSM2000'] = metaroot.style.magenta
    _color['ZprimeLH1000'] = metaroot.style.bright_green
    _color['ZprimeLH1250'] = metaroot.style.bright_red
    _color['ZprimeRH1000'] = metaroot.style.bright_green
    _color['ZprimeRH1250'] = metaroot.style.bright_red
    _color['ggH120']    = metaroot.style.bright_red
    _color['ggH125']    = metaroot.style.bright_red
    _color['ggH130']    = metaroot.style.bright_red
    _color['VBFH120']    = metaroot.style.bright_blue
    _color['VBFH125']    = metaroot.style.bright_blue
    _color['VBFH130']    = metaroot.style.bright_blue
    name = sample.name.split('.')[0]
    if 'Zprime' in name or 'ggH' in name or 'VBFH' in name:
        return metaroot.hist.PlotOptions(line_color = _color[name],
                                         line_width = 3,
                                         fill_color = _color[name],
                                         fill_style = 0 ) # hollow
    else:
        return metaroot.hist.PlotOptions(line_color = metaroot.style.black,
                                         fill_color = _color[name],
                                         fill_style = 1001 ) # solid

#_______________________________________________________________________________
def get_hists_and_options(hist_path, samples):
    hists = []
    plot_options = []
    labels = []
    for sample in samples:
        rename=sample.name
        h = sample.get(hist_path, rename=rename)
        if h:
            if ops.kfactors:
                k = get_kfactor(sample, hist_path)
                h.Scale(k)
            hists.append(h)
            plot_options.append(get_plot_option(sample))
            labels.append(sample.label)
        else:
            print "did not load %s from %s" % (sample.name, hist_path) #hack
    return hists, plot_options, labels

#_______________________________________________________________________________
def get_kfactor(sample, hist_path):
    k = 1.0
#    if sample.name.startswith('MultijetFakeTaus'):
#        k *= 2.0
#        
#    if sample.name.startswith('Zprime'):
#        k *= 1.4
    return k

#_______________________________________________________________________________
def is_left_sided(data_hist):
    data_mean = data_hist.GetMean()
    x_min = data_hist.GetXaxis().GetXmin()
    x_max = data_hist.GetXaxis().GetXmax()
    return data_mean < (x_max + x_min)/2.0

#______________________________________________________________________________
def calc_semi_data_driven_norm(data_hist, other_hists, name='', cache={}):
    if cache.has_key(name):
        return cache[name]
    else:
        N_data, N_data_error = metaroot.hist.calc_integral_and_error(data_hist)[0]
        N_other_and_errors = [ metaroot.hist.calc_integral_and_error(h) for h in other_hists ]
        N_other = sum( [ p[0] for p in N_other_and_errors ] )
#    k = kfactors.get(sample, hist_path)
#    if sample.name.startswith('Ztautau'):
#        k *= 1.3
#        
#    if sample.name.startswith('Zprime'):
#        k *= 1.4

    return k

#_______________________________________________________________________________
def is_left_sided(data_hist):
    data_mean = data_hist.GetMean()
    x_min = data_hist.GetXaxis().GetXmin()
    x_max = data_hist.GetXaxis().GetXmax()
    return data_mean < (x_max + x_min)/2.0

#______________________________________________________________________________
def calc_semi_data_driven_norm(data_hist, other_hists, name='', cache={}):
    if cache.has_key(name):
        return cache[name]
    else:
        N_data, N_data_error = metaroot.hist.calc_integral_and_error(data_hist)[0]
        N_other_and_errors = [ metaroot.hist.calc_integral_and_error(h) for h in other_hists ]
        N_other = sum( [ p[0] for p in N_other_and_errors ] )
        N_other_error = math.sqrt( sum( [ p[1]*p[1] for p in N_other_and_errors ] ) )
        N = N_data - N_other
        error = math.sqrt(N_data_error*N_data_error + N_other_error*N_other_error)
        cache[name] = N, error
    return N, error

#_______________________________________________________________________________
def do_the_great_rebinning(all_hists, hist_name, coarseness=1.0, log_x=False, log_y=False):
    ##-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    if hist_name in (
            'h_lep_pt',
            'h_el_pt',
            'h_mu_pt',
            'h_tau_pt',
            'h_mT',
            'h_mTtot',
            ):
        print 'rebinning ', hist_name
        for h in all_hists:
            metaroot.hist.rebin(h, bin_width=2.5*coarseness)
    ##-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    if hist_name in (
            'h_met_et',
            'h_vis_mass',
            'h_col_mass',
            'h_trans_eff_mass',
            'h_trans_mass',
            'h_eff_mass',
            'h_tau_lep_dpt',
            'h_mmc_mass',
            ):
        print 'rebinning ', hist_name
        for h in all_hists:
            metaroot.hist.rebin(h, bin_width=5.0*coarseness)
    ##-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    if hist_name in (
            'h_met_high',
            'h_mTtot_high',
            ):
        print 'rebinning ', hist_name
        for h in all_hists:
            metaroot.hist.rebin(h, bin_width=10*coarseness)
    ##-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    if hist_name in (
            'h_met_sumet',
            ):
        print 'rebinning ', hist_name
        for h in all_hists:
            metaroot.hist.rebin(h, bin_width=20*coarseness)
    ##-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    if hist_name in (
            'h_mjj',
            ):
        print 'rebinning ', hist_name
        for h in all_hists:
            metaroot.hist.rebin(h, bin_width=40*coarseness)
    ##-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    if hist_name in (
            'h_lep_eta',
            'h_el_eta',
            'h_mu_eta',
            'h_tau_eta',
            ):
        print 'rebinning ', hist_name
        for h in all_hists:
            metaroot.hist.rebin(h, bin_width=0.2*coarseness)
    ##-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    if hist_name in (
            'h_jj_deta',
            ):
        print 'rebinning ', hist_name
        for h in all_hists:
            metaroot.hist.rebin(h, bin_width=0.5*coarseness)
    ##-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    if hist_name in (
            'h_lep_phi',
            'h_el_phi',
            'h_mu_phi',
            'h_tau_phi',
            'h_met_phi',
            ):
        print 'rebinning ', hist_name
        for h in all_hists:
            metaroot.hist.rebin(h, n_merge=1*int(round(coarseness)))
    ##-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    if hist_name in (
            'h_taulep_dphi',
            'h_lepmet_dphi',
            'h_taumet_dphi',
            ):
        print 'rebinning ', hist_name
        for h in all_hists:
            metaroot.hist.rebin(h, n_merge=2*int(round(coarseness)))
    ##-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    return all_hists

#______________________________________________________________________________
def blind_the_data(data_hist, hist_name):
    if is_behind_blind(hist_name):
        print 'Blinding %s...' % hist_name
        h = data_hist.Clone(data_hist.GetName() + '_blinded')
        nbins = h.GetNbinsX()
        for i in xrange(nbins+2):
            x = h.GetBinLowEdge(i)
            if is_behind_blind(hist_name, x):
                h.SetBinContent(i, 0)
                h.SetBinError(i, 0.0)
        return h
    return data_hist

#______________________________________________________________________________
def is_behind_blind(hist_name, x=None):
    for key, val in blinding_dict.iteritems():
        if re.match(key, hist_name):
            if x is None or x >= val:
                return True
    return False

#______________________________________________________________________________
def draw_blind_line(hist_name, canv):
    canv.cd()
    stuff = {}
    for key,val in blinding_dict.iteritems():
        if re.match(key, hist_name):
            x = val
            break
    else:
         assert False
    stuff['blind_line'] = metaroot.utils.draw_vert_line(canv, x, color=ROOT.kWhite)
    stuff['blind_line_2'] = metaroot.utils.draw_vert_line(canv, x, color=ROOT.kBlack, style=7)
    lm = canv.GetLeftMargin()
    rm = 1. - canv.GetRightMargin()
    tm = 1. - canv.GetTopMargin()
    bm = canv.GetBottomMargin()
    xndc = (rm-lm)*((x-canv.GetUxmin())/(canv.GetUxmax()-canv.GetUxmin()))+lm
    blind_text = metaroot.utils.make_text(xndc+0.03, 0.60, 'Blinded', size=0.05, angle=90, color=ROOT.kBlack)
    blind_text.Draw()
#    blind_text_2 = metaroot.utils.make_text(xndc+0.03+0.002, 0.25-0.002, 'Blinded', size=0.05, angle=90, color=ROOT.kBlack)
#    blind_text_2.Draw()
    stuff['blind_text'] = blind_text
#    stuff['blind_text_2'] = blind_text_2
    return stuff

#______________________________________________________________________________
def set_hist_error_in_quadrature(h, error_hists):
    nbins = h.GetNbinsX()
    for i_bin in xrange(nbins+2):
        errors = [ he.GetBinError(i_bin) for he in error_hists ]
        h.SetBinError(i_bin, sum_in_quadrature(errors))

#______________________________________________________________________________
def set_hist_error_as_fraction_of_content(h, f):
    nbins = h.GetNbinsX()
    for i_bin in xrange(nbins+2):
        h.SetBinError(i_bin, h.GetBinContent(i_bin)*f)

#______________________________________________________________________________
def sum_in_quadrature(li):
    x = 0.0
    for y in li:
        x += y*y
    return math.sqrt(x)

#______________________________________________________________________________
def build_ratio_systs_hists(h_total, error_h_list):
    ratio_systs_hists = []
    nbins = h_total.GetNbinsX()
    for i_h, h in enumerate(error_h_list):
        h_ratio_band = h.Clone("h_ratio_band_%i" % i_h)
        set_hist_error_in_quadrature(h_ratio_band, error_h_list[:i_h+1])
        for i_bin in xrange(nbins+2):
            h_ratio_band.SetBinContent(i_bin, 1.0)
            c = h_total.GetBinContent(i_bin)
            e = h_ratio_band.GetBinError(i_bin) / c if c > 0.0 else 0.0
            h_ratio_band.SetBinError(i_bin, e)
        ratio_systs_hists.append(h_ratio_band)
    return ratio_systs_hists




#______________________________________________________________________________
if __name__ == '__main__': main()


