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
Follow the instruction given in the ***getting Start Guide*** in https://www.boost.org/. In sammary:

```bash
 cd $BUILD_DIR
 wget -c 'http://sourceforge.net/projects/boost/files/boost/1.72.0/boost_1_72_0.tar.bz2'
 tar -xvf boost_1_72_0.tar.bz2
 mv boost_1_72_0 boost
 cd boost
 ./bootstrap.sh --prefix=$(pwd)/boost --exec-prefix=$(pwd)/boost/bin --libdir=$(pwd)/boost/lib --includedir=$(pwd)/boost/include
 ./b2 install
```
Download, configure and install EIGEN
--------------------------------------

```bash
  cd $BUILD_DIR
  wget https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.tar.gz
  tar zxvf eigen-3.3.9.tar.gz
  mkdir eigen-3.3.9/build
  cd eigen-3.3.9/build
  cmake ../ -DCMAKE_INSTALL_PREFIX=build -DPYTHON_EXECUTABLE:FILEPATH=$BUILD_DIR/fenics2019_FSI/bin/python3
  make -j 8 install
```
Install MPI for python 
--------------------------------------
```bash
sudo apt-get install python3-pip
python3 -m pip install numpy mpi4py
```
Download, configure and install hdf5
-------------------------------------

```bash
  cd $BUILD_DIR
  cd boost
  wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.7/src/hdf5-1.10.7.tar.gz
  tar zxvf hdf5-1.10.7.tar.gz
  cd hdf5-1.10.7  
  ######### Note that this is one command split in several lines
  ./configure   --prefix=$BUILD_DIR/boost/hdf5-1.10.7_install \
  CC=mpicc CFLAGS=-O3 CXX=mpiCC CXXFLAGS=-O3 --enable-cxx \
  --enable-parallel --enable-unsupported
  
  make -j 8
  make install
  export LD_LIBRARY_PATH=${BUILD_DIR}/boost/hdf5-1.10.7_install/lib:$LD_LIBRARY_PATH
  export LD_RUN_PATH=${BUILD_DIR}/boost/hdf5-1.10.7_install/lib:$LD_RUN_PATH
 ```
Download, configure and install FEniCS python components
---------------------------------------------------------

```bash
  cd $BUILD_DIR
  wget https://bitbucket.org/fenics-project/ffc/downloads/ffc-2019.1.0.post0.tar.gz
  tar zxvf ffc-2019.1.0.post0.tar.gz
  cd ffc-2019.1.0.post0/
  sudo python3 setup.py install
```
