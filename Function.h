#ifndef Function_h
#define Function_h

#include <iostream>

class Function
{
	public:
		Function(const double & x_0, const int & n, const double & h);
		Function() = delete;
		~Function();
		
		double getX (const int i) const;
		double getF (const int i) const;
		double getX_0_ () const;
		int getN_ () const;
		double getH_ () const;
		void Sin ();
		
	
	private:
		double x_0_;
		int n_;
		double h_;
		double* x_i;
		double* f_i;
		
};

#endif