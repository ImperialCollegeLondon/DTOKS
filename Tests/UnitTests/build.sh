#!bin/bash

cd ../..
DTOKSUDIR=$PWD
UNITTESTDIR=$PWD/Tests/UnitTests/
cd $DTOKSUDIR

g++ -c -std=c++14 $DTOKSUDIR/Functions.cpp $DTOKSUDIR/threevector.cpp -I$DTOKSUDIR
cp Functions.o threevector.o $UNITTESTDIR/
cd $UNITTESTDIR

g++ -std=c++14 main.cpp Functions.o threevector.o -I$DTOKSUDIR -o main

rm Functions.o
