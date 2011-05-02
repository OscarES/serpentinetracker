#!/bin/bash

# preamble 
export BASE=`pwd`

# check OS (Serpentine tracker OS)
export STOS=`uname`

cd Dep 
tar xzvf xerces-c-3.0.1.tar.gz
cd xerces-c-3.0.1
export XERCESCROOT=`pwd`
./configure CC=gcc CXX=g++
make

cd ..
export UAP_FORTRAN_COMPILER=GFORTRAN  #(You may prefer to make this G95 instead.  Leave this variable unset to use ifort.) 
tar zxvf accelerator-ml.tar.gz
if [ "$STOS" == "Darwin" ]; then
    patch ./accelerator-ml/uap/trunk/Makefile ./patches/accelerator-ml-Makefile_macPatch
fi 
cd accelerator-ml/uap/trunk
export UAPROOT=`pwd`
make libs

cd $BASE/accformat/programs
export PYTHONVERSION=`echo \`python --version 2>&1\` | awk '{print $2}' | awk -F "." '{print $1 "." $2}'`
export PYTHONINCS=/usr/include/python$PYTHONVERSION
python setup.py build 

cd $BASE/Dep/
tar zxvf Minuit-1_7_9.tar.gz
cd $BASE/Dep/Minuit-1_7_9
export MINUITROOT=`pwd`
./configure
make 

cd $BASE/Dep
tar zxvf pyminuit-1.1.1.tgz
cd $BASE/Dep/pyminuit
Python setup.py build

cd $BASE/Dep
tar zxvf baseR3.14.10.tar.gz
cd $BASE/Dep/base-3.14.10
export EPICS_HOST_ARCH=`./startup/EpicsHostArch`
export EPICS_BASE=`pwd`
make
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$EPICS_BASE/lib/$EPICS_HOST_ARCH
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$EPICS_BASE/lib/$EPICS_HOST_ARCH
echo export EPICS_BASE=$EPICS_BASE
echo export EPICS_HOST_ARCH=$EPICS_HOST_ARCH
echo export PATH=$PATH:$EPICS_BASE/bin/$EPICS_HOST_ARCH/

cd $BASE/Dep/
tar zxvf pyepics-pyepics-3.1.2-0-g32c56ae.tar
cd $BASE/Dep/pyepics-pyepics-32c56ae/
python setup.py build


echo " "
echo "Please add the following lines to your shell configuration (.tcshrc .cshrc. .bashrc .profile)"
cd $BASE/accformat/programs
echo export PYTHONPATH=\$PYTHONPATH:$BASE/accformat/programs/build/`for i in \`ls build\`; do echo $i | grep lib; done`
cd $BASE/Dep/pyminuit
echo export PYTHONPATH=\$PYTHONPATH:$BASE/Dep/pyminuit/build/`for i in \`ls build\`; do echo $i | grep lib; done`
cd $BASE/Dep/pyepics-pyepics-32c56ae
echo export PYTHONPATH=\$PYTHONPATH:$BASE/Dep/pyepics-pyepics-32c56ae/build/`for i in \`ls build\`; do echo $i | grep lib; done`
echo export PYTHONPATH=\$PYTHONPATH:$BASE

echo export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$XERCESCROOT/src/.libs:$MINUITROOT/src/.libs:$EPICS_BASE/lib/$EPICS_HOST_ARCH
if [ "$STOS" == "Darwin" ]; then
    echo export DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:$XERCESCROOT/src/.libs:$MINUITROOT/src/.libs:$EPICS_BASE/lib/$EPICS_HOST_ARCH
fi 
echo export EPICS_BASE=$EPICS_BASE
echo export EPICS_HOST_ARCH=$EPICS_HOST_ARCH
echo export PATH=\$PATH:$EPICS_BASE/bin/$EPICS_HOST_ARCH/
