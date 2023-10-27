#!/bin/bash

update_catalogues=false 

while getopts ":u" opt; do
  case $opt in
    u)
      update_catalogues=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

echo "-Creating a conda environment"
conda deactivate
conda remove -y -n spear-mtb --all || :
conda env create -f ./environment.yml
source activate spear-mtb

echo "-Downloading the reference databases for CRyPTIC workflow"

if [ "$update_catalogues" = false ]; then
  wget -O="assets.tar.gz" https://figshare.com/ndownloader/files/38856204
  tar -xvzf assets.tar.gz
  # rm "assets.tar.gz"
fi

if [ "$update_catalogues" = true ]; then
  echo "-Updating the catalogues"
  cd ./assets/catalogues/NC_000962.3
  mv *.csv *.csv.BAK
  wget -O WHO_GARC1_v1.csv https://raw.githubusercontent.com/oxfordmmm/tuberculosis_amr_catalogues/public/catalogues/NC_000962.3/NC_000962.3_WHO-UCN-GTB-PCI-2021.7_v1.0_GARC1_RUS.csv
  wget -O CRyPTIC_GARC1_v1-311.csv https://raw.githubusercontent.com/oxfordmmm/tuberculosis_amr_catalogues/public/catalogues/NC_000962.3/NC_000962.3_CRyPTIC_v1.311_GARC1_RUS.csv
fi

echo "--Finished--"
