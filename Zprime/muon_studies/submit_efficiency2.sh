#!/bin/bash

cd /group/atlas/prj/leister/zpjun13/PennTau-ZprimeTauTauLepHad2012-00-00-00
source setup.sh -r
cd TestScripts
python muon_efficiency_subset2.py
