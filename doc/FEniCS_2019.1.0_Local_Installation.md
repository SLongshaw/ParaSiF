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

Download, configure and build pybind
-------------------------------------

```bash
  cd $BUILD_DIR
  export PYBIND11_VERSION=2.6.1
  wget https://github.com/pybind/pybind11/archive/v${PYBIND11_VERSION}.tar.gz
  tar zxvf v${PYBIND11_VERSION}.tar.gz
  mkdir pybind11-${PYBIND11_VERSION}/build
  cd pybind11-${PYBIND11_VERSION}/build
  cmake -DPYBIND11_TEST=off .. -DCMAKE_INSTALL_PREFIX=$(pwd)
  make install
```

Download, configure and prepare BOOST
--------------------------------------
Follow the instruction given in the ***getting Start Guide*** in https://www.boost.org/. In sammary, download boost file and uncobress in the required location. To keep all the installed libraries in one place, the recomnided location is $

```bash
 cd $BUILD_DIR
 wget -c 'http://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.bz2/download'
  tar -xvf boost_1_72_0.tar.bz2
mv boost_1_72_0 boost
cd boost
./bootstrap.sh --prefix=$(PWD)/boost 
./b2 install
```
