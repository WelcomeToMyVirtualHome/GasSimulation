#ifndef _GasSimulation_h_
#define _GasSimulation_h_

#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>
#include <stdlib.h>
#include <vector>
#include <limits>

#include "Utils.h"

namespace GasSimulation
{
	template<typename T>
	class Gas
	{
	public:
		Gas() 
		{ 
			srand48(time(NULL));
			F.resize((size_t)N*(N-1)/2); 
			r.reserve(N);
			p.resize(N);
			Utils::BaseVector<T> n_base(
				Utils::Vector<T>(a,0,0),
				Utils::Vector<T>(a/2,a*std::sqrt(3)/2,0),
				Utils::Vector<T>(a/2,a*std::sqrt(3)/6,a*std::sqrt(2.f/3)));
			base = n_base;	
		}
		Utils::BaseVector<T> base;
		size_t n = 5; // number of particles in one dimension
		size_t N = (int)pow(n,3); // number of all particles
		T T0 = 0; // [K]
		T a = 0.38; // [um]
		T L = 2.3; // [nm]
		T tau = 0.001; // dt [ps]
		const T m = 39.948; // [au]
		const T k_b = 0.00831;
		const T epsilon = 1;
		const T R = 0.38; // [nm]
		const T f = 10000;
		size_t s_d = 1000; // number of iterations of simulation
		const size_t s_0 = 1000; 
		const size_t s_xyz = 10; // xyz.dat: positions saved every 10 iterations
		const size_t s_out = 100; // out.dat: info saved every 100 iterations
		std::vector<Utils::Vector<T>> r;
		std::vector<Utils::Vector<T>> p;
		std::vector<Utils::Vector<T>> F;
		std::vector<Utils::Vector<T>> F_s;
		T V = 0;
		T T_temp = 0;
		T E_k = 0;
		T E_tot = 0;
		T P = 0;

		const char* PARAMS_FILE = "../Data/params.dat"; // parameters
		const char* R_FILE = "../Data/r.dat"; // initial positions
		const char* P_FILE = "../Data/p.dat"; // initial momentums
		const char* F_FILE = "../Data/F.dat"; // initial forces
		const char* XYZ_FILE = "../Data/xyz.dat"; // format "x y z", every chunk separated in file by 2 blank lines 
		const char* OUT_FILE = "../Data/out.dat"; // format "t V E_k E_tot T P"
		const char* EVAR_FILE = "../Data/E_var.dat"; // energy variance in function of tau
		const char* LACO_FILE = "../Data/lattice_constant.dat"; // potential energy in function of lattice constant
		const char* PT_FILE = "../Data/PT.dat";

		void SetParams(size_t n_n, T n_T0, T n_a)
		{
			n = n_n;
			T0 = n_T0;
			a = n_a;
			L = 1.22*n_a*(n_n-1);

			N = Utils::f_pow<size_t>(n,3);
			F.clear();
			F_s.clear();
			r.clear();
			p.clear();

			F.resize(N); 
			F_s.resize(N); 
			r.resize(N);
			p.resize(N);

			Utils::BaseVector<T> n_base(
				Utils::Vector<T>(a,0,0),
				Utils::Vector<T>(a/2,a*std::sqrt(3)/2,0),
				Utils::Vector<T>(a/2,a*std::sqrt(3)/6,a*std::sqrt(2.f/3)));
			base = n_base;	

			V = 0;
			T_temp = 0;
			E_k = 0;
			E_tot = 0;
			P = 0;
		}

		void CalculateInitialPositions()
		{	
			int index = 0;
			for(size_t i2 = 0; i2 < n; ++i2)
				for(size_t i1 = 0; i1 < n; ++i1)
					for(size_t i0 = 0; i0 < n; ++i0)
						r[index++] = base.b0*(i0 - (T)(n-1)/2) +base.b1*(i1 - (T)(n-1)/2) + base.b2*(i2 - (T)(n-1)/2); 
		}
		
		void CalculateInitialMomentum()	
		{
			Utils::Vector<T> p_sum;
			for(size_t i = 0; i < N; ++i)
			{
				Utils::Vector<T> p_i;
				p_i.x1 = std::sqrt(-m*k_b*T0*std::log(drand48()));
				p_i.x2 = std::sqrt(-m*k_b*T0*std::log(drand48()));
				p_i.x3 = std::sqrt(-m*k_b*T0*std::log(drand48()));		
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
			V = 0;
			P = 0;
			for(size_t i = 0; i < N; ++i){
				F[i].Reset();
			}

			for(size_t i = 0; i < N; ++i)
			{
				for(size_t j = 0; j < i && i > 0; ++j)
				{					
					T r_ij = r[i].Dist(r[j]);
					T V_p = epsilon*(Utils::f_pow<T>(R/r_ij,12)-2*Utils::f_pow<T>(R/r_ij,6));
					V += V_p;
					T multiplier = ((Utils::f_pow<T>(R/r_ij,12) - Utils::f_pow<T>(R/r_ij,6)) * 12*epsilon/Utils::f_pow<T>(r_ij,2));
					F[i] = F[i] + (r[i] + r[j]*(-1)) * multiplier;
					F[j] = F[j] + (r[j] + r[i]*(-1)) * multiplier;
				}
				T r_i = r[i].Norm();
				T V_s = r_i < L ? 0 : 0.5*f*Utils::f_pow<T>(r_i - L,2);	
				r_i < L ? F_s[i] = Utils::Vector<T>() : F_s[i] = r[i]*(f/r_i)*(L - r_i);
				V += V_s;
				F[i] = F[i] + F_s[i];
			}
		}

		void Simulation(T n_tau, T time)
		{
			CalculateInitialPositions();
			CalculateInitialMomentum(); 
			CalculatePotentialAndForces();
			std::ofstream outputXYZ;
			std::ofstream outputOUT;
			outputXYZ.open(XYZ_FILE, std::ios::trunc);
			outputOUT.open(OUT_FILE, std::ios::trunc);
			tau = n_tau;
			s_d = (size_t)(time/tau);
			for(size_t s = 0; s < s_d; ++s)
			{
				E_k = 0;
				for(size_t i = 0; i < N; ++i)
				{
					p[i] = p[i] + F[i]*(tau/2);
					r[i] = r[i] + p[i]*(tau/m);
					E_k += Utils::f_pow<T>(p[i].Norm(),2)/(2*m);
					P += F[i].Norm();
				}
				if(s % s_out == 0)
				{
					T_temp = 2/(3*N*k_b)*E_k;
					E_tot = E_k + V;
					P /= (4*M_PI*Utils::f_pow<T>(L,2));
					outputOUT << s*tau << " " << V  << " " << E_k << " " << E_tot << " " << T_temp << " " << P << "\n";
				}
				if(s % s_xyz == 0)
				{
					for(auto r_i : r)
						outputXYZ << r_i;
					outputXYZ << "\n\n";		
				}
				CalculatePotentialAndForces();
				for(size_t i = 0; i < N; ++i)
				{
					p[i] = p[i] + F[i]*(tau/2); 
				}
			}
			outputXYZ.close();	
			outputOUT.close();
		}

		T SimulationEnergyVariance(T tau, T time)
		{
			CalculateInitialPositions();
			CalculateInitialMomentum(); 
			CalculatePotentialAndForces();
			T E_total = 0;
			s_d = (size_t)(time/tau);
			std::vector<T> E_n(s_d);
			for(size_t s = 0; s < s_d; ++s)
			{
				E_k = 0;
				for(size_t i = 0; i < N; ++i)
				{
					p[i] = p[i] + F[i]*(tau/2);
					r[i] = r[i] + p[i]*(tau/m);
					E_k += pow(p[i].Norm(),2)/(2*m);
					P += F[i].Norm();
				}
				
				T_temp = 2/(3*N*k_b)*E_k;
				E_n[s] = E_k + V;
				E_total += E_n[s]; 
				CalculatePotentialAndForces();
				for(size_t i = 0; i < N; ++i)
				{
					p[i] = p[i] + F[i]*(tau/2); 
				}
			}
			T E_var = 0;
			E_total /= s_d;
			for(auto e : E_n)
				E_var = Utils::f_pow<T>(E_total-e,2);
			return E_var/s_d;
		}

		void StabilityTest(T tau1, T tau2, size_t numberOfPoints, T time)
		{
			std::ofstream output;
			output.open(EVAR_FILE, std::ios::trunc);
			T delta_tau = std::abs(tau1 - tau2)/numberOfPoints;
			tau = tau1;
			for(size_t i = 0; i < numberOfPoints; ++i)
			{
				output << tau << " " << SimulationEnergyVariance(tau,time) << "\n";
				tau += delta_tau;
			}
			output.close();
		}

		void CrystalConstantOptimization(T a1, T a2, size_t numberOfPoints)
		{
			std::ofstream output;
			output.open(LACO_FILE, std::ios::trunc);
			T delta_a = std::abs(a1 - a2)/numberOfPoints;
			a = a1;
			std::vector<T> V_n(numberOfPoints);
			for(size_t i = 0; i < numberOfPoints; ++i)
			{
				Utils::BaseVector<T> n_base(
					Utils::Vector<T>(a,0,0),
					Utils::Vector<T>(a/2,a*std::sqrt(3)/2,0),
					Utils::Vector<T>(a/2,a*std::sqrt(3)/6,a*std::sqrt(2.f/3)));
				base = n_base;	
				CalculateInitialPositions();
				CalculateInitialMomentum(); 
				CalculatePotentialAndForces();
				V_n[i] = V;
				output << a << " " << V << "\n";
				a += delta_a;
			}
			output.close();

			T V_min = std::numeric_limits<double>::max();
			size_t arg_min = 0;
			for(size_t i = 0; i < numberOfPoints; ++i)
				if(V_n[i] < V_min){
					V_min = V_n[i];
					arg_min = i;
				}
			T a_min = a1 + (arg_min+1)*delta_a;
			std::cout << "V_min= " << V_min << " a_min= " << a_min << '\n';
			a = a_min;
		}

		void Melting(T T1, T T2, T time, T n_tau)
		{
			T0 = T1;
			CalculateInitialPositions();
			CalculateInitialMomentum(); 
			CalculatePotentialAndForces();
			std::ofstream outputT;
			char buffer [50];
			sprintf(buffer, "./Data/melt_T0=%d.dat", (int)T0);
			outputT.open(std::string(buffer), std::ios::trunc);
			tau = n_tau;
			s_d = (size_t)(time/tau);
			T delta_T = (T1-T2)/s_d;
			for(size_t s = 0; s < s_d; ++s)
			{
				E_k = 0;
				for(size_t i = 0; i < N; ++i)
				{
					p[i] = p[i] + F[i]*(tau/2);
					r[i] = r[i] + p[i]*(tau/m);
					E_k += Utils::f_pow<T>(p[i].Norm(),2)/(2*m);
					P += F[i].Norm();
				}
				T0 += delta_T;
				T_temp = 2/(3*N*k_b)*E_k;
				E_tot = E_k + V;
				outputT << s*tau << " " << T_temp << "\n";
				CalculatePotentialAndForces();
				for(size_t i = 0; i < N; ++i)
				{
					p[i] = p[i] + F[i]*(tau/2); 
				}
			}
			outputT.close();	
		}

		void PressureTemperatureMonitor(T T1, T T2, size_t count, T time, T tau)
		{
			std::ofstream outputT, outputP, outputPT;
			outputPT.open(PT_FILE, std::ios::trunc);
			T delta_T = (T2-T1)/count;
			for(T t = T1; t < T2; t += delta_T)
			{
				T0 = t;
				CalculateInitialPositions();
				CalculateInitialMomentum(); 
				CalculatePotentialAndForces();
				char buffer [50];
				sprintf(buffer, "./Data/T_T=%d.dat", (int)T0);
				outputT.open(std::string(buffer), std::ios::trunc);
				sprintf(buffer, "./Data/P_T=%d.dat", (int)T0);
				outputP.open(std::string(buffer), std::ios::trunc);
				s_d = (size_t)(time/tau);
				T T_mean = 0;
				T P_mean = 0;
				for(size_t s = 0; s < s_d; ++s)
				{
					E_k = 0;
					for(size_t i = 0; i < N; ++i)
					{
						p[i] = p[i] + F[i]*(tau/2);
						r[i] = r[i] + p[i]*(tau/m);
						E_k += pow(p[i].Norm(),2)/(2*m);
						P += F_s[i].Norm();
					}
					T_temp = 2/(3*N*k_b)*E_k;
					P /= (4*M_PI*Utils::f_pow<T>(L,2));
					T t = s*tau;
					outputP << t << " " << P << "\n";
					outputT << t << " " << T_temp << "\n"; 
				
					T_mean += T_temp;
					P_mean += P;
					CalculatePotentialAndForces();
					for(size_t i = 0; i < N; ++i)
					{
						p[i] = p[i] + F[i]*(tau/2); 
					}
				}
				outputPT << P_mean/s_d << " " << T_mean/s_d << "\n";
				outputT.close();
				outputP.close();
			}
			outputPT.close();	
		}
		
		void FlushToFiles()
		{
			std::ofstream output;
			output.open(PARAMS_FILE, std::ios::trunc);
			output << n<<'\n'<<m<<'\n'<<epsilon<<'\n'<<R<<'\n'<<f<<'\n'<<L<<'\n'<<a<<'\n'<<T0<<'\n'<<tau<<'\n'<<s_0<<'\n'<<s_d<<'\n'<<s_out<<'\n'<<s_xyz;
			output.close();

			output.open(R_FILE, std::ios::trunc);
			for(auto r_i : r)
				output << r_i;
			output.close();

			output.open(P_FILE, std::ios::trunc);
			for(auto p_i : p)
				output << p_i;
			output.close();
			
			output.open(F_FILE, std::ios::trunc);
			for(auto F_i : F)
				output << F_i;
			output.close();
		}
	};
}

#endif
