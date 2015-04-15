#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
muon_resolution_plots.py

This is the method that writes and fills the plots for muon resolution studies
"""

import ROOT
import pyframe

GeV = 1000.0

def plot_mu_pt(alg, p, weight, prefix ='mu_', dir = ''):
    pTtruth = p.truth_pt
    reso = p.tlv.Pt()/pTtruth - 1.0
    alg.hist('g_%spT_resolution' % prefix, "ROOT.TH2F('$',';p_{T}(#mu)  [GeV];Resolution', 12, 0.0, 1200.0, 500, -1.5, 3.5)", dir=dir).Fill(pTtruth/GeV, reso)
