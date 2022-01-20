Install dependancies 
--------------------------------------------- 
Run the following to install **GNU compiler** and **make**
```bash
sudo apt update
sudo apt install build-essential
sudo apt install cmake
```
Download and build OpenMPI  as per https://www.open-mpi.org/software/ompi. After installation, it is also necessary to add the following line to the .bashrc file and source it :
```bash
export LD_LIBRARY_PATH=/usr/local/lib/
```


Create and set the installation and build folders
---------------------------------------------
```bash
  export INSTALL_FOLDER=`pwd`
  mkdir ${INSTALL_FOLDER}/FEniCS
  mkdir ${INSTALL_FOLDER}/FEniCS/V2019.1.0

  export BUILD_DIR=${INSTALL_FOLDER}/FEniCS/V2019.1.0
```

Load modules and set python paths/build paths
---------------------------------------------

```bash
  module load PrgEnv-gnu
  module load cray-python
  module load cmake

  cd $BUILD_DIR

  export PYTHONUSERBASE=${INSTALL_FOLDER}/.local
  export PATH=$PYTHONUSERBASE/bin:$PATH

  pip install --user virtualenv
  virtualenv --version
  virtualenv --system-site-packages fenics2019_FSI

  export PATH=$PATH:$BUILD_DIR/bin
  export PATH=$PATH:$BUILD_DIR/shared/bin
  export PYTHONPATH=$PYTHONPATH:$BUILD_DIR/lib/python3.8/site-packages
  export LD_LIBRARY_PATH=$BUILD_DIR/lib:$LD_LIBRARY_PATH
  export CC=cc
  export CXX=CC
```
