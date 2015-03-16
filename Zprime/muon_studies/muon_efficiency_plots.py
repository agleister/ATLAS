#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
muon_efficiency_plots.py

function for making muon pT plots for muon efficiency studies
"""
import ROOT

GeV = 1000.0

#_________________________________________________________________

def plot_true_mu_pt(alg, p, weight, prefix='true_mu_', dir = ''):
    alg.hist('h_%spT' % prefix, "ROOT.TH1F('$', ';p_{T}(#mu)  [GeV];Muons / (10 GeV)', 120, 0.0, 1200.0)", dir=dir).Fill(p.pt/GeV, weight)

def plot_mu_pt(alg, p, weight, prefix='mu_', dir = ''):
    alg.hist('h_%spT' % prefix, "ROOT.TH1F('$', ';p_{T}(#mu)  [GeV];Muons / (10 GeV)', 120, 0.0, 1200.0)", dir=dir).Fill(p.truth_pt/GeV, weight)
