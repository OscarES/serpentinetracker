#!/bin/bash

export BASE=`pwd`

cd Dep 
tar xzvf xerces-c-3.0.1.tar.gz
cd xerces-c-3.0.1
export XERCESCROOT=`pwd`
./configure CC=gcc CXX=g++
make

cd ..
export UAP_FORTRAN_COMPILER=GFORTRAN  #(You may prefer to make this G95 instead.  Leave this variable unset to use ifort.) 
tar zxvf accelerator-ml.tar.gz
cd accelerator-ml/uap/trunk
export UAPROOT=`pwd`
make libs

cd $BASE/accformat/programs
export PYTHONVERSION=`echo \`python --version 2>&1\` | awk '{print $2}' | awk -F "." '{print $1 "." $2}'`
export PYTHONINCS=/usr/include/python$PYTHONVERSION
python setup.py build 

echo " "
echo "Please add the following lines to your .bashrc"
echo export PYTHONPATH=\$PYTHONPATH:$BASE/accformat/programs/build/`for i in \`ls build\`; do echo $i | grep lib; done`
echo export PYTHONPATH=\$PYTHONPATH:$BASE
echo export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$XERCESCROOT/src/.libs

