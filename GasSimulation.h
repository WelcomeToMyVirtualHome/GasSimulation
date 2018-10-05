#ifndef _GasSimulation_h_
#define _GasSimulation_h_

#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <vector>

#include "Utils.h"

namespace GasSimulation
{
	template<typename T>
	class Gas
	{
	public:
		Gas() { srand48(time(NULL)); }
		Utils::BaseVector<T> base;
		const int n = 5;
		const double T0 = 100; // [K]
		const double m = 39.948; // [au]
		const double kb = 0.00831;
		const double a = 0.38; // [um]
		
		void CalculateInitialPositions(std::string outputFileName)
		{	
			std::ofstream output;
			output.open(outputFileName);
			Utils::Vector<T> r;
			for(int i2 = 0; i2 < n; ++i2)
			{
				for(int i1 = 0; i1 < n; ++i1)
				{
					for(int i0 = 0; i0 < n; ++i0)
					{
						Utils::Vector<T> r0 = base.b0*(double)(i0 - (n-1)/2); 
						Utils::Vector<T> r1 = base.b1*(double)(i1 - (n-1)/2); 
						Utils::Vector<T> r2 = base.b2*(double)(i2 - (n-1)/2); 
						r = r0 + r1 + r2;
						output << r;
					}
				}	
			}
			output.close();	
		}
		
		void CalculateInitialMomentum(std::string outputFileName)	
		{
			std::vector<Utils::Vector<T>> momentums;
			const int N = (int)pow(n,3);
			momentums.reserve(N);
			Utils::Vector<T> momentumSum;
			for(int i = 0; i < N; ++i)
			{
				Utils::Vector<T> p;
				p.x1 = sqrt(-m*kb*T0*log(drand48()));
				p.x2 = sqrt(-m*kb*T0*log(drand48()));
				p.x3 = sqrt(-m*kb*T0*log(drand48()));		
				p.x1 *= drand48() > 1./2  ? -1 : +1;
				p.x2 *= drand48() > 1./2  ? -1 : +1;
				p.x3 *= drand48() > 1./2  ? -1 : +1;
				momentums.push_back(p);	
				momentumSum = momentumSum + p;
			}
			std::cout << "P= " << momentumSum << "p_dash= " << sqrt(m*kb*T0) << "\n";
			std::ofstream output;
			output.open(outputFileName);
			for(auto &p : momentums)
			{
				p = p  + momentumSum*(-1./N);
				output << p;
			}
			output.close();	
		}
	};
}

#endif