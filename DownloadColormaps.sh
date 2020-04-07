#!/bin/sh

curl -O http://www.fabiocrameri.ch/resources/lapaz.zip
curl -O http://www.fabiocrameri.ch/resources/vik.zip
curl -O http://www.fabiocrameri.ch/resources/turku.zip

unzip lapaz.zip
unzip vik.zip
unzip turku.zip

rm -r lapaz.zip
rm -r vik.zip
rm -r turku.zip

echo "If you publish a figure with the colormap, make sure to cite the following papers"
echo "Crameri, F. (2018). Scientific colour-maps. Zenodo. http://doi.org/10.5281/zenodo.1243862"

echo "Crameri, F. (2018), Geodynamic diagnostics, scientific visualisation and StagLab 3.0, Geosci. Model Dev., 11, 2541-2562, doi:10.5194/gmd-11-2541-2018"

