#include <iostream>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;


// costanti numeriche
constexpr double L = 500.0;
constexpr double x0 = 200.0;
constexpr double q = 2.0;
constexpr double sigma = 20.0;
constexpr double a = 250.0;
constexpr double b = 260.0;
constexpr double V0 = 1.7;
constexpr int Nx = 2000;
constexpr double dt = 0.1;
constexpr int nTimeSteps = 2000;
constexpr double hx = L/Nx;



// funzione che risolve un'equazione di matrice tridiagonale
void tridiagonal_solve(int n, complex<double>* upper, complex<double>* lower, complex<double>* diag, 
								complex<double>* sol, complex<double>* b);



int main() {
	
	// numero di iterazioni eseguite
	int N_it = 0;
	
	
	// output file
	fstream fout;
	fout.open("TDSchroedinger.dat",ios::out);
	fout.precision(10);
	
	
	// funzioni d'onda nuova (n2) e vecchia (n1) COMPLESSE
	complex<double>* psi_n1 = new complex<double>[Nx];
	complex<double>* psi_n2 = new complex<double>[Nx];
	complex<double>* psi_tmp = new complex<double>[Nx-2];
	
	
	// inizializzo psi_0
	psi_n2[0] = complex<double>(0.0,0.0);
	psi_n2[Nx-1] = complex<double>(0.0,0.0);
	for(int j=1; j<=Nx-2; j++)
	{
		// definita come funzione complessa
		psi_n2[j] = polar(1.0, q*j*hx) * exp(-0.5*pow((j*hx-x0)/sigma, 2.0));
		
	}
	
	
	// calcolo norma psi_0 con il metodo di Simpson
	double Norm_psi0 = 0.0;
	
	Norm_psi0 += (1.0/3.0)*hx*norm(psi_n2[0]);
	for(int j=1; j<=Nx-2; j++)
	{
		if(j%2==0) Norm_psi0 += (2.0/3.0)*hx*norm(psi_n2[j]);
		else Norm_psi0 += (4.0/3.0)*hx*norm(psi_n2[j]);
		
	}
	Norm_psi0 += (1.0/3.0)*hx*norm(psi_n2[Nx-1]);
	Norm_psi0 = sqrt(Norm_psi0);
	
	
	// normalizzo psi_0
	for(int j=0; j<Nx; j++)
	{
		psi_n2[j] /= Norm_psi0;
		
	}
	
	
	// vettori operativi per calcolo matriciale
	// contenenti solo le componenti non banali
	complex<double>* up = new complex<double>[Nx-2];
	complex<double>* low = new complex<double>[Nx-2];
	complex<double>* diag = new complex<double>[Nx-2];
	complex<double>* b_k = new complex<double>[Nx-2];
	
	double Vj;
	
	
	
	// calolo ricorsivo
	while(N_it < nTimeSteps)
	{
		N_it++;
		
		// sposto funzione vecchia
		for(int j=0; j<Nx; j++)
		{
			psi_n1[j] = psi_n2[j];
		
		}
		
		
		// -- preparo vettori per la funzione matriciale --
		
		// up & low
		up[0] = complex<double>(1.0,0.0);
		up[Nx-3] = complex<double>(0.0,0.0);
		low[0] = complex<double>(0.0,0.0);
		low[Nx-3] = complex<double>(1.0,0.0);
		for(int j=1; j<=Nx-4; j++)
		{
			up[j] = complex<double>(1.0,0.0);
			low[j] = complex<double>(1.0,0.0);
		}
		
		// diag
		diag[0] = complex<double>(-2.0-2.0*hx*hx*0, 4.0*hx*hx/dt);
		diag[Nx-3] = complex<double>(-2.0-2.0*hx*hx*0, 4.0*hx*hx/dt);
		for(int j=1; j<=Nx-4; j++)
		{
			// calcolo Vj
			if(a <= (j+1)*hx && (j+1)*hx <= b) Vj = V0;
			else Vj = 0.0;
			
			diag[j] = complex<double>(-2.0-2.0*hx*hx*Vj, 4.0*hx*hx/dt);
		}
		
		// b_k
		for(int j=0; j<=Nx-3; j++)
		{
			// calcolo Vj
			if(a <= (j+1)*hx && (j+1)*hx <= b) Vj = V0;
			else Vj = 0.0;
			
			b_k[j] = complex<double>(2.0+2.0*hx*hx*Vj, 4.0*hx*hx/dt) * psi_n1[j+1] - psi_n1[(j+1)-1] - psi_n1[(j+1)+1];
		}
		
		
		// risolvo matrice
		tridiagonal_solve(Nx-2, up, low, diag, psi_tmp, b_k);
		
		
		
		// sistemo risultati nel vettore finale correttamente
		psi_n2[0] = complex<double>(0.0,0.0);
		psi_n2[Nx-1] = complex<double>(0.0,0.0);
		for(int j=0; j<=Nx-3; j++) psi_n2[j+1] = psi_tmp[j];
		
		
		// check norma
		if(N_it%200 == 0)
		{
			// calcolo norma psi_0 con il metodo di Simpson
			Norm_psi0 = 0.0;
	
			Norm_psi0 += (1.0/3.0)*hx*norm(psi_n2[0]);
			for(int j=1; j<=Nx-2; j++)
			{
				if(j%2==0) Norm_psi0 += (2.0/3.0)*hx*norm(psi_n2[j]);
				else Norm_psi0 += (4.0/3.0)*hx*norm(psi_n2[j]);
		
			}
			Norm_psi0 += (1.0/3.0)*hx*norm(psi_n2[Nx-1]);
			Norm_psi0 = sqrt(Norm_psi0);
			
			
			cout << "Ripetizione " << N_it << ": norma = " << Norm_psi0 << endl;
			
		}
		
		
		
		// scrivi su file
        if(N_it % 10 == 0) {
        	
			for(int j=0; j<Nx; j++)
				fout << hx*j << " " << abs(psi_n2[j]) << endl;
			
			fout << endl << endl;
			
        }
   
		
	}
	
	
	fout.close();
	return 0;
	
	
	
}