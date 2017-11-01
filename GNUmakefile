#------------------------------------------------------------------------------

#SOURCE= main.cpp Breakup.cpp Constants.cpp DTOKSU.cpp threevector.cpp Functions.cpp Models/Model.cpp Models/ChargingModel.cpp Models/HeatingModel.cpp Models/ForceModel.cpp  Models/PlasmaGrid.cpp Elements/Matter.cpp Elements/Iron.cpp Elements/Beryllium.cpp Elements/Tungsten.cpp Elements/Graphite.cpp
SOURCE= *.cpp Models/*.cpp Elements/*.cpp
INC=-I. -I./Models -I./Elements
CC=g++
CFLAGS= -std=c++14

#------------------------------------------------------------------------------



all: main

main: $(SOURCE)
	$(CC) $(CFLAGS) $(INC) -o $@ $^

.PHONY: clean

clean:
	rm main
