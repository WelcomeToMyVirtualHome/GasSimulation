#include <iostream>
#include "GasSimulation.h"
#include <cmath>

int main()
{
	GasSimulation::Gas<double> gas;
	Utils::BaseVector<double> base(Utils::Vector<double>(gas.a,0,0),Utils::Vector<double>(gas.a/2,gas.a*sqrt(3)/2,0),Utils::Vector<double>(gas.a/2,gas.a*sqrt(3)/6,gas.a*sqrt(2.f/3)));
	gas.base = base;	
	gas.CalculateInitialPositions();
	gas.CalculateInitialMomentum(); 
	gas.CalculatePotentialAndForces();
	gas.FlushToFiles();
	gas.Simulation();
	return 1;
}
