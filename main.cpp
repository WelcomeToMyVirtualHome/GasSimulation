#include <iostream>
#include "GasSimulation.h"
#include <cmath>

int main()
{
	GasSimulation::Gas<double> gas;
	gas.SetParams(5, 100., 0.38, 2.3);
	gas.CalculateInitialPositions();
	gas.CalculateInitialMomentum(); 
	gas.CalculatePotentialAndForces();
	gas.FlushToFiles();	
	// gas.Simulation();

	// gas.SetParams(3, 1000., 0.38, 1.2);
	// gas.StabilityTest(std::pow(10,-6), 5*std::pow(10,-3), 100, 1.);
	
	// gas.SetParams(5, 1000., 0.38, 2.3);
	// gas.CrystalConstantTest(0.35, 0.40, 1000);



	return 1;
}
