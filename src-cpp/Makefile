all: shoal

Optimizer.o: Optimizer.cpp 
	g++ -I include -std=c++11 -O3 -c Optimizer.cpp
asa103.o: asa103.cpp 
	g++ -I include -std=c++11 -O3 -c asa103.cpp
shoal.o: shoal.cpp
	g++ -I include -std=c++11 -O3 -c shoal.cpp 
shoal: shoal.o asa103.o Optimizer.o
	g++ -std=c++11 -O3 shoal.o asa103.o Optimizer.o -o shoal -lpthread
