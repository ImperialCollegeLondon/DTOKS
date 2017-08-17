#!bin/bash

cd ../../..
DTOKSUDIR=$PWD
ELEMENTSDIR=$DTOKSUDIR/Elements
MODELSDIR=$DTOKSUDIR/Models
INTEGRATIONTESTDIR=$DTOKSUDIR/Tests/IntegrationTests/Forcing
cd $DTOKSUDIR

g++ -c -std=c++14 $INTEGRATIONTESTDIR/main.cpp MathHeader.cpp DTOKSU.cpp Constants.cpp threevector.cpp Functions.cpp Models/Model.cpp Models/ForceModel.cpp  Models/PlasmaGrid.cpp Models/Model.cpp Elements/Matter.cpp Elements/Iron.cpp Elements/Beryllium.cpp Elements/Tungsten.cpp Elements/Graphite.cpp -I. -I./Models -I./Elements

mv *.o $INTEGRATIONTESTDIR

cd $INTEGRATIONTESTDIR

g++ -std=c++14 main.cpp MathHeader.o Constants.o threevector.o Model.o PlasmaGrid.o ForceModel.o Beryllium.o Iron.o Tungsten.o Graphite.o Matter.o Functions.o -o main -I$DTOKSUDIR -I$INTEGRATIONTESTDIR -I$MODELSDIR -I$ELEMENTSDIR

rm main.o Beryllium.o Iron.o Tungsten.o Graphite.o Matter.o Functions.o Model.o PlasmaGrid.o threevector.o ForceModel.o Constants.o MathHeader.o
