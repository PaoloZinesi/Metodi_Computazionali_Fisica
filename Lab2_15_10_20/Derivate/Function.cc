#include <iostream>
#include "math.h"

#include "Function.h"
	
// se non viene specificato nulla la memoria per le derivate NON viene allocata
Function::Function(const double & x_0, const int & n, const double & h):
	x_0_(x_0),
	n_(n),
	h_(h),
	alloc_f1_(false) {
		
		x_i = new double[n_];
		f_i = new double[n_];
		
		
		// inizializza griglia
		for (int k=0; k<n_; k++)
		{
			x_i[k] = x_0_ + k * h_;
		}
		
		
	}
	

// costruttore per allocare memoria per le derivate	
Function::Function(const double & x_0, const int & n, const double & h, const bool alloc):
	Function(x_0, n, h) {
	
	if(alloc)
	{
		alloc_f1_ = true;
		
		
		// alloca elementi della matrice
		f1_i = new double*[4];
		for(int k=0; k<4; k++)
		{
			f1_i[k] = new double[n_];
		}
		
		
	}

}	
	

Function::~Function() {
			
		delete[] x_i;
		delete[] f_i;
		
		
		// elimina solo se la memoria è stata allocata
		if(alloc_f1_)
		{
			for(int k=0; k<4; k++)
			{
				delete[] f1_i[k];
			}
			delete[] f1_i;
			
		}
			
	}
	
	
double Function::getX (const int i) const
{
	return x_i[i];
}


double Function::getF (const int i) const
{
	return f_i[i];
}


double Function::getF1 (const int metodo, const int i) const
{
	return f1_i[metodo][i];
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

void Function::Cos ()
{
	for(int i=0; i<n_; i++)
	{
		f_i[i] = cos(x_i[i]);
	}
}


// calcola la derivata per il metodo scelto
void Function::Derivative (const int tipo)
{
	switch (tipo) {
		
	  case 0: 		// derivata destra O(h)
	    
	  	for(int i=0; i<(n_-1); i++)
	  	{
		  	f1_i[0][i] = (f_i[i+1]-f_i[i])/h_;
		
	  	}
		
		// calcolo anche i valori marginali
		f1_i[0][n_] = f1_i[0][n_-1]; 
		
	  	break;
	  
	  
  	  case 1: 		// derivata sinistra O(h)
	    
  	  	for(int i=1; i<n_; i++)
  	  	{
  		  	f1_i[1][i] = (f_i[i]-f_i[i-1])/h_;
		
  	  	}
		
  		// calcolo anche i valori marginali
  		f1_i[1][0] = f1_i[1][1]; 
		
  	  	break;
	  
	  
  	  case 2: 		// derivata O(h^2)
	    
  	  	for(int i=1; i<(n_-1); i++)
  	  	{
  		  	f1_i[2][i] = (f_i[i+1] - f_i[i-1]) / (2 * h_);
		
  	  	}
		
  		// calcolo anche i valori marginali
  		f1_i[2][0] = (f_i[1]-f_i[0])/h_;
		f1_i[2][n_] = (f_i[n_]-f_i[n_-1])/h_;
		
  	  	break;
	  
	  
  	  case 3: 		// derivata O(h^4)
	    
  	  	for(int i=2; i<(n_-2); i++)
  	  	{
  		  	f1_i[3][i] = (-f_i[i+2] + 8*f_i[i+1] - 8*f_i[i-1] + f_i[i-2]) / (12 * h_);
		
  	  	}
		
  		// calcolo anche i valori marginali
  		f1_i[3][0] = (f_i[1]-f_i[0])/h_;
		f1_i[3][1] = (f_i[2] - f_i[0]) / (2 * h_);
		f1_i[3][n_-1] = (f_i[n_] - f_i[n_-2]) / (2 * h_);
		f1_i[3][n_] = (f_i[n_]-f_i[n_-1])/h_;
		
  	  	break;
		
		
	}
	
	return;
	
}

