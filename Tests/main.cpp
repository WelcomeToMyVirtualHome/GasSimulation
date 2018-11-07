#include <iostream>
#include "../libs/GasSimulation.h"

int main(int argc, char **argv)
{
	GasSimulation::Gas<double> gas;
	gas.SetParams(5, 1000., 0.38);
	gas.CalculateInitialPositions();
	gas.CalculateInitialMomentum(); 
	gas.CalculatePotentialAndForces();
	gas.FlushToFiles();	
	
	gas.SetParams(3, 1000., 0.38);
	gas.StabilityTest(std::pow(10,-6), 5*std::pow(10,-3), 100, 1.);
	
	gas.SetParams(5, 2000., 0.38);
	gas.CrystalConstantOptimization(0.35, 0.40, 1000);

	gas.PressureTemperatureMonitor(500, 2500, 4, 5., 0.001);
	return 1;
}
