# thmdyn

Scripts for producing data and figures related to:

Anheuser, J., Liu, Y., and Key, J.: A climatology of thermodynamic vs. dynamic Arctic wintertime sea ice thickness effects during the CryoSat-2 era, submitted to: The Cryosphere. 2022

Data can be found at:



Scripts utilize Python 3.8.3 and Jupyter. Python packages required include:

numpy
pandas
glob
datetime
scipy
xarray
pyproj
matplotlib
cartopy
cmocean
h5py
pyhdf

Auxilliary data required include:

AMSR-E and AMSR2 brightness temperatures: https://doi.org/10.5067/AMSR-E/AE_SI25.003; https://doi.org/10.5067/TRUIAL3WPAUP
AMSR-E and AMSR2 SIC: https://doi.org/10.5067/AMSR-E/AE_SI25.003; https://doi.org/10.5067/TRUIAL3WPAUP
AWI CS2SMOS: https://www.meereisportal.de
Sea ice motion vectors: https://doi.org/10.5067/INAWUWO7QH7B
MOSAiC drift track:  https://doi.pangaea.de/10.1594/PANGAEA.937193
