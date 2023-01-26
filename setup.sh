#!/usr/bin/env bash

echo "1.Creating a conda environment"
conda deactivate
conda remove -y -n spear-mtb --all || :
conda env create -f ./environment.yml
source activate spear-mtb

echo "2.Downloading the reference databases for clockwork workflow"

wget -O="assets.tar.gz" https://figshare.com/ndownloader/files/38856204
tar -xvzf assets.tar.gz
# rm "assets.tar.gz"

echo "3.Updating the catalogues"
cd ./assets/catalogues/NC_000962.3
mv *.csv *.csv.BAK
wget -O WHO-GARC1.csv https://raw.githubusercontent.com/oxfordmmm/tuberculosis_amr_catalogues/public/catalogues/NC_000962.3/NC_000962.3_WHO-UCN-GTB-PCI-2021.7_v1.0_GARC1_RUS.csv