#!/bin/bash
set -e

./blatHuman.bash
R CMD BATCH removeHuman.R
./blastNt.bash
R CMD BATCH analyzeNt.R
