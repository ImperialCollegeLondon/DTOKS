#!bin/bash

cd ../..
DTOKSUDIR=$PWD
CHARGINGDIR=$PWD/Tests/UnitTests/Charging
UNITTESTDIR=$PWD/Tests/UnitTests/
BOOSTDIR=$PWD/boost_1_64_0/
cd $DTOKSUDIR

g++ -c -std=c++14 $DTOKSUDIR/Constants.cpp $DTOKSUDIR/Functions.cpp $DTOKSUDIR/threevector.cpp -I$DTOKSUDIR
mv Functions.o Constants.o threevector.o $UNITTESTDIR/
cd $UNITTESTDIR

g++ -std=c++14 main.cpp Constants.o Functions.o threevector.o -I$DTOKSUDIR -I$CHARGINGDIR -I$BOOSTDIR -o main

rm Functions.o threevector.o Constants.o
