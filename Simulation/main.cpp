#include <iostream>
#include "../libs/GasSimulation.h"

int main(int argc, char **argv)
{
	if (argc < 4)
	{
		std::cout<<"usage: n_atoms_axis T_0 lattice_constant dt time"<<'\n';	
		return 0;
	}
	GasSimulation::Gas<double> gas;
	gas.SetParams(atof(argv[1]), atof(argv[2]), atof(argv[3]));
	gas.Simulation(atof(argv[4]), atof(argv[5]));
	return 1;
}
