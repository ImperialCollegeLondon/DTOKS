#!bin/bash
 
g++ -std=c++14 main.cpp DTOKSU.cpp ChargingEquation/ChargingModel.cpp HeatingEquation/HeatingModel.cpp ForceEquation/ForceModel.cpp Elements/Matter.cpp -o main -I. -I./HeatingEquation -I./ForceEquation -I./ChargingEquation -I./Elements
