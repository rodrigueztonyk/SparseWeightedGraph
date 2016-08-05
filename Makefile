# This is a Makefile for using nauty things

GPATH	= /opt/gurobi651/linux64/
NPATH	= /opt/nauty26r4/

CC 	= gcc
CPP	= g++

CPPARG	= -m64 -g -lpthread -lm -O3

GI	= -I$(GPATH)include/
GL	= -L$(GPATH)lib/ -lgurobi_c++ -lgurobi65
NI	= -I$(NPATH) $(NPATH)nauty.a

SWG 	= SparseWeightedGraph.cpp $(NI) -w 
SS	= SparseWeightedGraph.cpp SparseWeightedGraph.hpp

dumb : dumb.cpp $(SS)
	$(CPP) --std=c++11 $(SWG) dumb.cpp -o dumb -w
