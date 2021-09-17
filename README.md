# ParMupSiF - Parallel Partitioned Multi-physical Simulation Framework

Parallel Partitioned Multi-physical Simulation Framework is developed based on the MUI library. It offers a platform where users can carry out multi-physical (mainly fluid-structure interaction) studies using supercomputers.

The framework uses a partitioned approach to couple two or more physical domains together for multi-physics simulations. It takes several advantages of the MUI library:

• Flexibility to select different solvers for each physical domain;

• Flexibility to extend the number of physical domains;

• Good scalability on communications among physical domains for large simulations;

• Keeping the development of the coupled solvers decoupled. It allows for easier independent testing of each solver and avoids potential incompatibilities between two solvers (i.e if they both use a certain library X each one depends on a different version of it which are not compatible);

• Coupling of multiple solvers which have different programming language interfaces (e.g C, C++, FORTRAN and Python);

• "Plug and play" strategy. One solver can be replaced by another incrementally and without the need of recompiling if the MUI interface is used as a common adaptor;

• Use of multiple solvers which have two incompatible licenses exploiting the dual licensing of the MUI library (both solvers are never mixed source-wise or binary-wise).

**This framework is a beta software at the moment and under active development** to involve more solvers as well as more physical domains. Such infrastructure will make it possible to simulate large multi-physical problems and simulate complicated multi-physical cases by using supercomputing facilities.

## Licensing

Copyright (C) 2021 Engineering and Environment Group, Scientific Computing Department, Science and Technology Facilities Council, UK Research and Innovation. All rights reserved.

This code is licensed under the GNU General Public License version 3

The Parallel Partitioned Multi-physical Simulation Framework provides FEniCS v2019.1.0 <-> MUI v1.1.2 <-> OpenFOAM v6 two-way coupling.

## Acknowledgements
The Parallel Partitioned Multi-physical Simulation Framework is developed at the [Scientific Computing Department](https://www.scd.stfc.ac.uk/) of the [Science and Technology Facilities Council](https://stfc.ukri.org/). If you this framework, please cite us:

*Liu, W., Wang, W., Skillen, A., Longshaw, S.M., Moulinec, C. and Emerson, D.R. (2021). A Parallel Partitioned Approach on Fluid-Structure Interaction Simulation Using the Multiscale Universal Interface Coupling Library. In: 14th World Congress In Computational Mechanics (WCCM) And ECCOMAS Congress 2020.*

## Publication

Liu, W., Longshaw, S.M., Skillen, A., Emerson, D.R., Valente, C. and Gambioli, F. (2021). A High-Performance Open-Source Solution for Multiphase Fluid-Structure Interaction. In: 31st International Ocean and Polar Engineering Conference (ISOPE).

Liu, W., Wang, W., Skillen, A., Longshaw, S.M., Moulinec, C. and Emerson, D.R. (2021). A Parallel Partitioned Approach on Fluid-Structure Interaction Simulation Using the Multiscale Universal Interface Coupling Library. In: 14th World Congress In Computational Mechanics (WCCM) And ECCOMAS Congress 2020.

Liu, W., Wang, W., Skillen, A., Fernandez, E.R., Longshaw, S. and Sawko, R. (2020). Code Development on Parallel Partitioned Fluid-Structure Interaction Simulations. Project Report. [online] STFC e-Pub. Available at: https://epubs.stfc.ac.uk/work/46400815.

## Contact

Should you have any question please do not hesitate to contact the developers

## Installation

**Step One: Install FEniCS v2019.0.1**

Please follow the instruction from FEniCS Prjoect Manual to install the FEniCS version 2019.1.0, which was released on April 19th 2019: https://fenicsproject.org/download/

If building FEniCS from source, please refer to the **Stable version** section to install the latest stable release of FEniCS (version 2019.1.0) from [this link](https://fenics.readthedocs.io/en/latest/installation.html#from-source)

**Step Two: Obtain MUI and install FSI coupling lab and wrappers**

```bash
cd parMupSiF/src/MUI_Utility
wget https://github.com/MxUI/MUI/archive/1.1.1.tar.gz
tar -xzf 1.1.1.tar.gz
cd MUI-1.1.1/
cp -r ../couplingFSILab ./
cd ..
rm 1.1.1.tar.gz
mv MUI-1.1.1/ ../../../MUI

cd ../../../MUI/wrappers/C
make

cd ../Python
```
In here, we need to change the Makefile in the MUI/wrappers/Python folder, by relpace Line 5 from "USE_RBF = False" to "USE_RBF = True"
Then:

```
make COMPILER=GCC package
make pip-install

cd ../../couplingFSILab/wrappers/C
make -f Makefile_CAPI

cd ../Python
make COMPILER=GCC package
make pip-install
cd ../../../../parMupSiF/src/MUI_Utility
```

For more information on MUI, please refer to the MxUI GitHub organisation pages https://github.com/MxUI

**Step Three: Install OpenFOAM v6 with MUI patches**

```bash
cd ../../../
mkdir -p OpenFOAM/v6
cd OpenFOAM/v6

wget -O - http://dl.openfoam.org/source/6 | tar xvz
wget -O - http://dl.openfoam.org/third-party/6 | tar xvz

mv OpenFOAM-6-version-6 OpenFOAM-6
mv ThirdParty-6-version-6 ThirdParty-6

cp ../../parMupSiF/src/MUI_Utility/OpenFOAM_patch/* ./

./patch_OF6-MUI
```

Please place a symlink with absolute path to the MUI library into the new folder "ThirdParty-6/MUI" before move forward.

Please note: 1. make sure the symlink to the MUI library is not relative path as a script takes a copy of these files; 2. Whenever MUI has been updated, the old platform files of OpenFOAM need to be cleaned and re-run Allwmake to update the copy that gets used by the OpenFOAM applications.

For more information, please refere to parMupSiF/src/MUI_Utility/OpenFOAM_patch/README, OpenFOAM-6/README-MUI and ThirdParty-6/README-MUI

```bash
source OpenFOAM-6/etc/bashrc
cd OpenFOAM-6


./Allwmake

cd ../../../parMupSiF/src/MUI_Utility
```
For more information on OpenFOAM v6, please refer to the OpenFOAM Foundation: https://openfoam.org/download/6-source/

**Step Four: Install parMupSiF OpenFOAM solvers**

```bash
cd ../CFD/OpenFOAM/V6/applications/solvers/

# Source OpenFOAM run functions
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# Compile parMupSiF solvers
wmake pimpleFSIFoam
wmake interFSIFoam
```

## Source and export before run parMupSiF cases

```bash
source /path/to/dolfin/dolfin.conf
source /path/to/OpenFOAM-6/etc/bashrc
export PYTHONPATH= /path/to/parMupSiF/src/CSM/FEniCS/V2019.1.0:$PYTHONPATH
```

## Demos

Demos of both the single phase FSI case (3-D flow past a flexible beam) and the multiphase FSI case (3-D dam break with a flexible beam) are in parMupSiF/demo/FSIBeam_pimpleFSIFoam and parMupSiF/demo/FSIBeam_interFSIFoam, respectively.
