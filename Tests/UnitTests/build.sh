#!bin/bash

cd ../..
DTOKSUDIR=$PWD
CHARGINGDIR=$PWD/Tests/UnitTests/Charging
HEATINGDIR=$PWD/Tests/UnitTests/Heating
UNITTESTDIR=$PWD/Tests/UnitTests/
BOOSTDIR=$PWD/boost_1_64_0/
cd $DTOKSUDIR

g++ -c -std=c++14 $DTOKSUDIR/Constants.cpp $DTOKSUDIR/MathHeader.cpp $DTOKSUDIR/Functions.cpp $DTOKSUDIR/threevector.cpp -I$DTOKSUDIR
mv Functions.o Constants.o MathHeader.o threevector.o $UNITTESTDIR/
cd $UNITTESTDIR

g++ -std=c++14 main.cpp Constants.o MathHeader.o Functions.o threevector.o -I$DTOKSUDIR -I$CHARGINGDIR -I$HEATINGDIR -I$BOOSTDIR -o main

rm Functions.o threevector.o Constants.o MathHeader.o
