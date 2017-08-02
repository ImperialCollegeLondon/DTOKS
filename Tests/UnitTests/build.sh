#!bin/bash

cd ../..
DTOKSUDIR=$PWD
CHARGINGDIR=$PWD/Tests/UnitTests/Charging
UNITTESTDIR=$PWD/Tests/UnitTests/
BOOSTDIR=/home/ls5115/Software/boost_1_64_0/
cd $DTOKSUDIR

g++ -c -std=c++14 $DTOKSUDIR/Functions.cpp $DTOKSUDIR/threevector.cpp -I$DTOKSUDIR
mv Functions.o threevector.o $UNITTESTDIR/
cd $UNITTESTDIR

g++ -std=c++14 main.cpp Functions.o threevector.o -I$DTOKSUDIR -I$CHARGINGDIR -I$BOOSTDIR -o main

rm Functions.o
