#!bin/bash

cd ../..
DTOKSUDIR=$PWD
ERRORESTIMATEDIR=$PWD/Tests/ModelTesting/
cd $DTOKSUDIR

g++ -c -std=c++14 threevector.cpp Constants.cpp MathHeader.cpp Models/PlasmaGrid.cpp Models/Model.cpp Models/ForceModel.cpp Models/HeatingModel.cpp Elements/Beryllium.cpp Elements/Iron.cpp Elements/Tungsten.cpp Elements/Graphite.cpp Elements/Matter.cpp Functions.cpp -I$DTOKSUDIR -I$DTOKSUDIR/Elements -I$DTOKSUDIR/Models


mv *.o $ERRORESTIMATEDIR

cd $ERRORESTIMATEDIR

g++ -std=c++14 main.cpp threevector.o Constants.o Model.o MathHeader.o ForceModel.o PlasmaGrid.o HeatingModel.o Beryllium.o Iron.o Tungsten.o Graphite.o Matter.o Functions.o -o main -I$DTOKSUDIR -I$ERRORESTIMATEDIR -I$DTOKSUDIR/Models -I$DTOKSUDIR/Elements

rm HeatingModel.o Beryllium.o Iron.o Tungsten.o Graphite.o Matter.o Functions.o Constants.o threevector.o Model.o PlasmaGrid.o ForceModel.o MathHeader.o
