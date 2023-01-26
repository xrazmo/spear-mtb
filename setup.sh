#!/usr/bin/env bash

conda deactivate
conda remove -y -n spear-mtb --all || :
conda env create -f ./environment.yml
source activate spear-mtb
