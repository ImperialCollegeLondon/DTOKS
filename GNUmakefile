#------------------------------------------------------------------------------

SOURCE= *.cpp Models/*.cpp Elements/*.cpp
INC=-I. -I./Models -I./Elements
LINKNETCDF=-L${NETCDFDIR}/lib -lnetcdf_c++ -lhdf5_hl -lhdf5 -lz -lm
CC=g++
CFLAGS= -std=c++14

#------------------------------------------------------------------------------



all: main

main: $(SOURCE)
	$(CC) $(CFLAGS) $(INC) -o $@ $^ -L/usr/local/lib -lnetcdf_c++ -lhdf5_hl -lhdf5 -lz -lm

.PHONY: clean

clean:
	rm main
