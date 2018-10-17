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
		const int n = 5; // number of particles in one dimension
		const int N = (int)pow(n,3); // number of all particles
		const T T0 = 100; // [K]
		const T m = 39.948; // [au]
		const T k_b = 0.00831;
		const T a = 0.38; // [um]
		const T epsilon = 1;
		const T R = 0.38; // [nm]
		const T L = 2.3; // [nm]
		const T f = 10000;
		const T tau = 0.01; // dt [ps]
		const int s_d = 1000; // number of iterations of simulation
		const int s_0 = 1000; 
		const int s_xyz = 10; // xyz.dat: positions saved every 10 iterations
		const int s_out = 100; // out.dat: info saved every 100 iterations
		std::vector<Utils::Vector<T>> r;
		std::vector<Utils::Vector<T>> p;
		std::vector<Utils::Vector<T>> F;
		T V;
		T T_temp;
		T E_k;
		T E_tot;
		T P;
		T t;

		const char* R_FILE = "r.dat"; // initial positions
		const char* P_FILE = "p.dat"; // initial momentums
		const char* F_FILE = "F.dat"; // initial forces
		const char* XYZ_FILE = "xyz.dat"; // format "x y z", every chunk separated in file by 2 blank lines 
		const char* OUT_FILE = "out.dat"; // format "t V E_k E_tot T P"

		void CalculateInitialPositions()
		{	
			Utils::Vector<T> r_i;
			r.reserve(N);
			for(int i2 = 0; i2 < n; ++i2)
			{
				for(int i1 = 0; i1 < n; ++i1)
				{
					for(int i0 = 0; i0 < n; ++i0)
					{
						Utils::Vector<T> r0 = base.b0*(double)(i0 - (n-1)/2); 
						Utils::Vector<T> r1 = base.b1*(double)(i1 - (n-1)/2); 
						Utils::Vector<T> r2 = base.b2*(double)(i2 - (n-1)/2); 
						r_i = r0 + r1 + r2;
						r.push_back(r_i);
					}
				}	
			}
		}
		
		void CalculateInitialMomentum()	
		{
			p.resize(N);
			Utils::Vector<T> p_sum;
			for(int i = 0; i < N; ++i)
			{
				Utils::Vector<T> p_i;
				p_i.x1 = sqrt(-m*k_b*T0*log(drand48()));
				p_i.x2 = sqrt(-m*k_b*T0*log(drand48()));
				p_i.x3 = sqrt(-m*k_b*T0*log(drand48()));		
				p_i.x1 *= drand48() > 1./2  ? -1 : +1;
				p_i.x2 *= drand48() > 1./2  ? -1 : +1;
				p_i.x3 *= drand48() > 1./2  ? -1 : +1;	
				p_sum = p_sum + p_i;
				p[i] = p_i;
			}
			for(auto &p_i : p)
			{
				p_i = p_i + p_sum*(-1./N);
			}
		}

		void CalculatePotentialAndForces()
		{
			F.resize((size_t)N*(N-1)/2);
			V = 0;
			for(int i = 0; i < N; ++i)
			{
				T r_i = r[i].Norm();
				T V_s = r_i < L ? 0 : 0.5*f*pow(r_i - L,2);	
				r_i < L ? F[i] = Utils::Vector<T>() : F[i] = F[i] + r[i]*(f/r_i)*(L - r_i);
				V += V_s;
				for(int j = 0; j < i && i > 0; ++j)
				{					
					T r_ij = r[i].Dist(r[j]);
					T V_p = epsilon*(pow(R/r_ij,12)-2*pow(R/r_ij,6));
					V += V_p;
					T multiplier = ((pow(R/r_ij,12) - pow(R/r_ij,6)) * 12*epsilon*pow(r_ij,-2));
					F[i] = F[i] + (r[i] + r[j]*(-1)) * multiplier;
					F[j] = F[j] + (r[j] + r[i]*(-1)) * multiplier;
				}
			}
		}

		void Simulation()
		{
			std::ofstream outputXYZ;
			std::ofstream outputOUT;
			outputXYZ.open(XYZ_FILE, std::ios::trunc);
			outputOUT.open(OUT_FILE, std::ios::trunc);
			t = 0;
			for(int s = 0; s < s_d; ++s)
			{
				E_k = 0;
				P = 0;
				for(int i = 0; i < N; ++i)
				{
					p[i] = p[i] + F[i]*(tau/2);
					r[i] = r[i] + p[i]*(tau/m);
					E_k += pow(p[i].Norm(),2)/(2*m);
					P += F[i].Norm();
				}
				if(s % s_out == 0)
				{
					T_temp = 2/(3*N*k_b)*E_k;
					E_tot = E_k + V;
					t += s*tau;
					P /= (4*M_PI*pow(L,2));
					outputOUT << t << " " << V  << " " << E_k << " " << E_tot << " " << T_temp << " " << P << "\n";
				}
				if(s % s_xyz == 0)
				{
					for(auto r_i : r)
						outputXYZ << r_i;
					outputXYZ << "\n\n";		
				}
				CalculatePotentialAndForces();
				for(int i = 0; i < N; ++i)
				{
					p[i] = p[i] + F[i]*(tau/2); 
				}

			}
			outputXYZ.close();	
			outputOUT.close();
		}

		
		void FlushToFiles()
		{
			std::ofstream output;
			output.open(R_FILE);
			for(auto r_i : r)
				output << r_i;
			output.close();

			output.open(P_FILE);
			for(auto p_i : p)
				output << p_i;
			output.close();
			
			output.open(F_FILE);
			for(auto F_i : F)
				output << F_i;
			output.close();
		}
	};
}

#endif
