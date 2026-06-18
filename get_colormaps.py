import os
import shutil
import urllib.request
import zipfile

url = "https://zenodo.org/records/5501399/files/ScientificColourMaps7.zip?download=1"
zip_file = "ScientificColourMaps7.zip"
wanted = ["lapaz", "vik", "turku"]

print("Downloading ScientificColourMaps7...")
urllib.request.urlretrieve(url, zip_file)

print("Unzipping...")
with zipfile.ZipFile(zip_file, "r") as z:
    z.extractall()

for name in wanted:
    src = os.path.join("ScientificColourMaps7", name)
    dst = name

    if os.path.exists(dst):
        shutil.rmtree(dst)

    shutil.move(src, dst)
    print(f"Moved {name}")

shutil.rmtree("ScientificColourMaps7")

print()
print("Done!")
print()
print("If you publish a figure with the colormap, make sure to cite:")
print("Crameri, F. (2018). Scientific colour-maps. Zenodo. doi:10.5281/zenodo.1243862")
print("Crameri, F. (2018), Geodynamic diagnostics, scientific visualisation and StagLab 3.0, Geosci. Model Dev., 11, 2541-2562.")
