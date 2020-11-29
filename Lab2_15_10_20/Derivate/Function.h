#ifndef Function_h
#define Function_h

#include <iostream>

class Function
{
	public:
		Function(const double & x_0, const int & n, const double & h);
		Function(const double & x_0, const int & n, const double & h, const bool alloc);
		Function() = delete;
		~Function();
		
		double getX (const int i) const;
		double getF (const int i) const;
		double getF1 (const int metodo, const int i) const;
		
		double getX_0_ () const;
		int getN_ () const;
		double getH_ () const;
		
		
		void Sin ();
		void Cos ();
		void Derivative (const int tipo);
		
		
	
	private:
		double x_0_;
		int n_;
		double h_;
		double* x_i;		// griglia di punti
		double* f_i;		// valori della funzione
		
		bool alloc_f1_;		// matrice allocata true/false
		double** f1_i;		// derivata prima della funzione nei 4 metodi rispettivi
		
};

#endif