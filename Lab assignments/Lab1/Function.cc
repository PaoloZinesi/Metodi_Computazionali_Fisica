#include <iostream>
#include "math.h"

#include "Function.h"
	

Function::Function(const double & x_0, const int & n, const double & h):
	x_0_(x_0),
	n_(n),
	h_(h) {
		
		x_i = new double[n_];
		f_i = new double[n_];
		
		for (int k=0; k<n_; k++)
		{
			x_i[k] = x_0_ + k * h_;
		}
		
	}
	

Function::~Function(){
			
		delete[] x_i;
		delete[] f_i;
			
	}
	
	
double Function::getX (const int i) const
{
	return x_i[i];
}


double Function::getF (const int i) const
{
	return f_i[i];
}


double Function::getX_0_ () const
{
	return x_0_;
}

int Function::getN_ () const
{
	return n_;
}

double Function::getH_ () const
{
	return h_;
}

void Function::Sin ()
{
	for(int i=0; i<n_; i++)
	{
		f_i[i] = sin(x_i[i]);
	}
}

