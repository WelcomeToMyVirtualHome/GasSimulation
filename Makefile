CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic

GasSimulation:	main.o GasSimulation.o
		$(CXX) -o GasSimulation main.o GasSimulation.o $(CXXFLAGS) $(LIBS)
main.o:	main.cpp
		$(CXX) -o main.o -c main.cpp $(CXXFLAGS)
GasSimulation.o:	GasSimulation.cpp
		$(CXX) -o GasSimulation.o -c GasSimulation.cpp $(CXXFLAGS)
clean:
		rm *.o
