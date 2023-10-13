# Stimulated Brillouin Spectroscopy (SBS)

This library is designed to process data acquired from Stimulated Brillouin Spectroscopy (SBS) experiments.
It provides a user-friendly interface for analyzing SBS data and includes various functions for fitting, visualization and processing of data. 
The use of Matlab ensures compatibility with a wide range of data formats and processing tools, while the incorporation of gpufit allows for fast and efficient fitting of complex models even for large datasets. 
The processing is GPU-accelerated (based on Gpufit), so a compatible Nvidia GPU-card (and CUDA) is needed to run this software.
The library is designed to be flexible and easily extendable, making it a valuable tool for researchers in the field of SBS.

This repository contains data and code used by Fan Yang during his work on the paper. Most of the code was written by Sebastian Hambura.

## Repository organisation 
### Matlab-pipeline
This is the software used to process and analyze the recorded data. Read more about it [here](./Matlab-pipeline/readme.md).
To try it, go to *SI-scripts > douple-peak-zebrafish* : it contains the dataset and processing scripts of 2 datasets (XY and XZ plane of a zebrafish).

The code and the precompiled binaries were tested and are running with the following configuration : Matlab 2021.b and Cuda 11.5. If you're running a different setup, you might have to recompile or rewrite a bit the code.

### Gpufit-addition
This folder contains all the Gpufit customization needed to process SBS data. [Read this file](./gpufit-addition/readme.md) for more information about how to build the custom Gpufit version. You can also [download this fork](https://github.com/prevedel-lab/Gpufit), already containing the customizations. 

### SI-scripts
This regroups different scripts and dataset used in the supplementary figures. The processing scripts are sometimes a bit different compared to the normal processing. [Click here](./SI-scripts/readme.md) for more information.

## References
- Przybylski, A., Thiel, B., Keller-Findeisen, J., Stock, B., and Bates, M. (2017). Gpufit: An open -source toolkit for GPU-accelerated curve fitting. Sci. Rep. 7, 15722
- Yoon-Oh Tak (2023). Multipage TIFF stack (https://www.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack), MATLAB Central File Exchange. Retrieved October 13, 2023.  

## Contact
- [Prevedel Lab, EMBL](https://www.prevedel.embl.de/)
- [Link to Fan](https://people.ucas.edu.cn/~yangfan)

## Licence
Copyright (C) 2022-2023 Sebastian Hambura

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.