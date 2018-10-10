#ifndef _Utils_h_
#define _Utils_h_

#include <iostream>
#include <fstream>

namespace Utils
{
	template <typename T>
	class Vector
	{
	public:
		Vector() : x1(0), x2(0), x3(0) {}
		Vector(T a, T b, T c) : x1(a), x2(b), x3(c) {} 
		T x1, x2, x3;
		
		Vector operator +(const Vector &v){ return Vector(x1 + v.x1, x2 + v.x2, x3 + v.x3); }
		Vector operator *(double a){ return Vector(a*x1, a*x2, a*x3); }
		T Dist(const Vector &v){ return sqrt(pow(x1 - v.x1,2) + pow(x2 - v.x2,2) + pow(x3 - v.x3,2) ); }
		T Norm() { return sqrt(pow(x1,2) + pow(x2,2) + pow(x3,2)); }
		void Reset() { x1 = 0; x2 = 0; x3 = 0; }
		friend std::ostream &operator <<(std::ostream &s, const Vector &v ){ return s  << v.x1 << " " << v.x2 << " " << v.x3 << "\n"; }
	};

	template <typename T>
	class BaseVector
	{
	public:	
		BaseVector() {}
		BaseVector(Vector<T> a, Vector<T> b, Vector<T> c) : b0(a), b1(b), b2(c) {}
		Vector<T> b0,b1,b2;
		friend std::ostream &operator <<(std::ostream &s, const BaseVector &b) { return s << b.b0 << b.b1 << b.b2; }
	};

	template <typename T>
	static std::vector<Utils::Vector<T>> ReadFileVectors(std::string filename)
	{
		std::ifstream input;
		input.open(filename);
		std::vector<Utils::Vector<T>> vectors;	
		T val1, val2, val3;
		while(input>>val1 && input>>val2 && input>>val3)
		{
			vectors.push_back(Utils::Vector<T>(val1,val2,val3));
		}
		input.close();
		return vectors;
	}		
}

#endif
