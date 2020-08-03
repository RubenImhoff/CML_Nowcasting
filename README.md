# Rainfall Nowcasting Using Commercial Microwave Links and Radar Rainfall QPE
This repository contains the Python and R scripts that were used in the publication "Rainfall Nowcasting Using Commercial Microwave Links" by Imhoff et al. (2020a). The repository contains five main folders: 'Adjusted_pySTEPS_scripts', 'Data', 'GridPointsToGeotiff', 'Nowcasting' and 'RAINLINK_Scripts'.

The Commercial Microwave Link (CML) data is first interpolated to a 1 km2 grid that is the same as the grid of the KNMI radar rainfall data (https://doi.org/10.4121/uuid:05a7abc4-8f74-43f4-b8b1-7ed7f5629a01). RAINLINK is employed for this purpose (Overeem et al., 2013, 2016). The used R-scripts of RAINLINK are present in the folder 'RAINLINK_Scripts' and all necessary (CML) data is provided in the folder 'Data' (i.e., the parameter values, the interpolation grid and the output CML data). After this interpolation step, the data is available as 15-minute rainfall intensity per grid cell location. Subsequently, the Python script in 'GridPointsToGeotiff' is used to save the data in a GeoTIFF-format that has the same grid and extent as the KNMI radar rainfall data. 

The radar and CML rainfall estimates are available on the same grid now. Probabilistic nowcasts are made with pysteps (Pulkkinen et al., 2019) with the scripts in the folder 'Nowcasting'. These scripts follow the method described in Imhoff et al. (2020b). The folder 'Adjusted_pySTEPS_scripts' contains the scripts that were adjusted, as compared to the original model configurations of pysteps (Pulkkinen et al., 2019), to allow for nowcasting with the GeoTIFF files. For the distribution of the latest pysteps version, see also https://github.com/pySTEPS/pysteps. 

# Used algorithms
RAINLINK (Overeem et al., 2013, 2016) is available via: https://github.com/overeem11/RAINLINK.
pysteps (Pulkkinen et al., 2019) is available via: https://github.com/pySTEPS/pysteps. 

# Usage
For the use of pysteps, please see the aforementioned github page. Example run scripts are also available on the github page of pysteps. 

# Terms of use
pysteps (Copyright (c) 2019, PySteps developers) is provided under the BSD 3-Clause and RAINLINK (Copyright (C) 2019 Aart Overeem) under a GNU General Public License. Both algorithms are free to use. The pysteps and RAINLINK folders in this repository follow the mentioned terms of use for both algorithms. This repository is in no way a new means of distributing the source code of both algorithms. It is meant for the reproducibility of the study conducted by Imhoff et al. (2020a). For the latest distribution and examples of both algorithms, please see the aforementioned github pages. 

The scripts and data in the other folders ('Data', 'GridPointsToGeotiff', 'Nowcasting') are free to use: you can redistribute them and/or modify them under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your opinion) any later version. These scripts are redistributed in the hope that they will be useful, but without any warranty and without even the implied warranty of merchantability or fitness for a particular purpose. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licences/. The producer of this tool is by no means to be held responsible for the use of and the results of this tool.

# References
Imhoff, R. O., Overeem, A., Brauer, C. C., Leijnse, H., Weerts, A. H., & Uijlenhoet, R. (2020a). Rainfall Nowcasting Using Commercial Microwave Links. Manuscript submitted to Geophysical Research Letters.

Imhoff, R. O., Brauer, C. C., Overeem, A., Weerts, A. H., & Uijlenhoet, R. (2020b). Spatial and temporal evaluation of radar rainfall nowcasting techniques on 1,533 events. Water Resources Research, 56, e2019WR026723. doi: 10.1029/2019WR026723.

Overeem, A., Leijnse, H., & Uijlenhoet, R. (2013). Country-wide rainfall maps from cellular communication networks. Proceedings of the National Academy of Sciences, 110 (8), 2741-2745. doi: 10.1073/pnas.1217961110.

Overeem, A., Leijnse, H., & Uijlenhoet, R. (2016a). Retrieval algorithm for rainfall mapping from microwave links in a cellular communication network. Atmospheric Measurement Techniques, 9 (5), 2425-2444. doi: 10.5194/amt-9-2425-2016.

Pulkkinen, S., Nerini, D., Pérez Hortal, A. A., Velasco-Forero, C., Seed, A., Germann, U., & Foresti, L. (2019). Pysteps: an open-source Python library for probabilistic precipitation nowcasting (v1.0). Geoscientific Model Development, 12 (10), 4185-4219. doi: 10.5194/gmd-12-4185-2019. 

