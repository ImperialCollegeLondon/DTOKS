#!bin/bash

cd ../../..
DTOKSUDIR=$PWD
ELEMENTSDIR=$DTOKSUDIR/Elements
MODELSDIR=$DTOKSUDIR/Models
INTEGRATIONTESTDIR=$DTOKSUDIR/Tests/IntegrationTests/Charging/
cd $DTOKSUDIR

g++ -c -std=c++14 $INTEGRATIONTESTDIR/main.cpp Constants.cpp threevector.cpp Functions.cpp Models/Model.cpp Models/ChargingModel.cpp  Models/PlasmaGrid.cpp Models/Model.cpp Elements/Matter.cpp Elements/Iron.cpp Elements/Beryllium.cpp Elements/Tungsten.cpp Elements/Graphite.cpp -I. -I./Models -I./Elements

mv *.o $INTEGRATIONTESTDIR

cd $INTEGRATIONTESTDIR

g++ -std=c++14 main.cpp Constants.o threevector.o Model.o PlasmaGrid.o ChargingModel.o Beryllium.o Iron.o Tungsten.o Graphite.o Matter.o Functions.o -o main -I$DTOKSUDIR -I$INTEGRATIONTESTDIR -I$MODELSDIR -I$ELEMENTSDIR

rm main.o Beryllium.o Iron.o Constants.o Tungsten.o Graphite.o Matter.o Functions.o Model.o PlasmaGrid.o threevector.o ChargingModel.o
