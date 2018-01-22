#------------------------------------------------------------------------------

SOURCE= *.cpp Models/*.cpp Elements/*.cpp
INC=-I. -I./Models -I./Elements -I/home/ls5115/Software/config4cpp/include/
CC=g++
CFLAGS= -std=c++14

#------------------------------------------------------------------------------


all: main

main: $(SOURCE)
	$(CC) $(CFLAGS) $(INC) -o $@ $^ -L/usr/local/lib -L/home/ls5115/Software/config4cpp/lib -lnetcdf_c++ -lhdf5_hl -lhdf5 -lz -lm -lconfig4cpp
#-lconfig4cpp.a


.PHONY: clean

clean:
	rm main
