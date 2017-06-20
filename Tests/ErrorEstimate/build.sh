#!bin/bash

HEATINGMATTERDIR=$HOME/Software/HeatingMatter/
ERRORESTIMATEDIR=$HOME/Software/HeatingMatter/Tests/ErrorEstimate/
cd $HEATINGMATTERDIR

g++ -c -std=c++14 main.cpp HeatingModel.cpp Beryllium.cpp Iron.cpp Tungsten.cpp Graphite.cpp Matter.cpp Emissiv.cpp RefractiveIndexGen.cpp Functions.cpp

gfortran -c MIEV0.f MVTstOld.f ErrPack.f

mv *.o $ERRORESTIMATEDIR

cd $ERRORESTIMATEDIR

g++ -std=c++14 main.cpp HeatingModel.o Beryllium.o Iron.o Tungsten.o Graphite.o Matter.o Emissiv.o Functions.o RefractiveIndexGen.o ErrPack.o MIEV0.o MVTstOld.o -o main -lgfortran -I$HEATINGMATTERDIR -I$ERRORESTIMATEDIR

rm main.o HeatingModel.o Beryllium.o Iron.o Tungsten.o Graphite.o Matter.o Emissiv.o Functions.o RefractiveIndexGen.o ErrPack.o MIEV0.o MVTstOld.o
