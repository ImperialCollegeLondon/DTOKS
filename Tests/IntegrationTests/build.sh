#!bin/bash

DTOKSUDIR=/home/ls5115/Software/DTOKS-U/
ELEMENTSDIR=/home/ls5115/Software/DTOKS-U/Elements
MODELSDIR=/home/ls5115/Software/DTOKS-U/Models
INTEGRATIONTESTDIR=/home/ls5115/Software/DTOKS-U/Tests/IntegrationTests/
cd $DTOKSUDIR

g++ -c -std=c++14 $INTEGRATIONTESTDIR/main.cpp DTOKSU.cpp threevector.cpp Functions.cpp Models/Model.cpp Models/ChargingModel.cpp Models/HeatingModel.cpp Models/ForceModel.cpp  Models/PlasmaGrid.cpp Models/Model.cpp Elements/Matter.cpp Elements/Iron.cpp Elements/Beryllium.cpp Elements/Tungsten.cpp Elements/Graphite.cpp -I. -I./Models -I./Elements

mv *.o $INTEGRATIONTESTDIR

cd $INTEGRATIONTESTDIR

g++ -std=c++14 main.cpp threevector.o Model.o PlasmaGrid.o HeatingModel.o Beryllium.o Iron.o Tungsten.o Graphite.o Matter.o Functions.o -o main -I$DTOKSUDIR -I$INTEGRATIONTESTDIR -I$MODELSDIR -I$ELEMENTSDIR

rm main.o HeatingModel.o Beryllium.o Iron.o Tungsten.o Graphite.o Matter.o Functions.o Model.o PlasmaGrid.o
