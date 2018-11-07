CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic -O3

GasSimulation:	main.o
		$(CXX) -o GasSimulation main.o $(CXXFLAGS) $(LIBS)
main.o:	main.cpp
		$(CXX) -o main.o -c main.cpp $(CXXFLAGS)
clean:
		rm *.o
