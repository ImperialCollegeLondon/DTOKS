#!bin/bash
 
g++ -std=c++14 main.cpp DTOKSU.cpp threevector.cpp Functions.cpp Models/Model.cpp Models/ChargingModel.cpp Models/HeatingModel.cpp Models/ForceModel.cpp Elements/Matter.cpp Elements/Iron.cpp Elements/Beryllium.cpp Elements/Tungsten.cpp Elements/Graphite.cpp -o main -I. -I./Models -I./Elements
