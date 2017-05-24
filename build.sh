#!bin/bash
 
g++ -std=c++14 main.cpp DTOKSU.cpp threevector.cpp Functions.cpp ChargingEquation/ChargingModel.cpp HeatingEquation/HeatingModel.cpp ForceEquation/ForceModel.cpp Elements/Matter.cpp Elements/Iron.cpp Elements/Beryllium.cpp Elements/Tungsten.cpp Elements/Graphite.cpp -o main -I. -I./HeatingEquation -I./ForceEquation -I./ChargingEquation -I./Elements
