#include <iostream>
#include "GasSimulation.h"
#include <cmath>

int main()
{
	GasSimulation::Gas<double> gas;
	double a = gas.a;
	Utils::Vector<double> b0(a,0,0);
	Utils::Vector<double> b1(a/2,a*sqrt(3)/2,0);
	Utils::Vector<double> b2(a/2,a*sqrt(3)/6,a*sqrt(2.f/3));
	Utils::BaseVector<double> base(b0,b1,b2);
	gas.base = base;	

	gas.CalculateInitialPositions();
	gas.CalculateInitialMomentum(); 
	gas.CalculatePotentialAndForces();
	gas.FlushToFiles();
	return 1;
}
