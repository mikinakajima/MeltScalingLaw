#!/bin/sh


https://zenodo.org/records/5501399#.YlSGbzfMIqs
unzip ScientificColourMaps7.zip 

unzip lapaz.zip
unzip vik.zip
unzip turku.zip

mv ScientificColourMaps7/lapaz .
mv ScientificColourMaps7/vik .
mv ScientificColourMaps7/turku .

rm -r ScientificColourMaps7


echo "If you publish a figure with the colormap, make sure to cite the following papers"
echo "Crameri, F. (2018). Scientific colour-maps. Zenodo. http://doi.org/10.5281/zenodo.1243862"
echo "Crameri, F. (2018), Geodynamic diagnostics, scientific visualisation and StagLab 3.0, Geosci. Model Dev., 11, 2541-2562, doi:10.5194/gmd-11-2541-2018"

