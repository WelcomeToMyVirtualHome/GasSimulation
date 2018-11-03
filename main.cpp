#include <iostream>
#include "GasSimulation.h"
#include <cmath>

int main()
{
	GasSimulation::Gas<double> gas;
	gas.SetParams(5, 1000., 0.38, 2.3);
	gas.CalculateInitialPositions();
	gas.CalculateInitialMomentum(); 
	gas.CalculatePotentialAndForces();
	gas.FlushToFiles();
	
	gas.Simulation();
	
	return 1;
}
