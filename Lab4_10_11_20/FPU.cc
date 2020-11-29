#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;


int main(int argc, char* argv[])
{
	
	int N, N_modi, N_it_max, N_cycl_max, N_it=0;
	double m = 1.0, k = 1.0;
	double alpha, beta;
	double ampl = 1.0;
	double dt;
	
	cout << "Inserisci N (punti materiali) ";
	cin >> N;
	cout << "Inserisci N_modi (modi normali da considerare) ";
	cin >> N_modi;
	cout << "Inserisci N_cycl_max (numero massimo di cicli rispetto alla frequenza fondamentale) ";
	cin >> N_cycl_max;
	cout << "Inserisci alpha ";
	cin >> alpha;
	cout << "Inserisci beta ";
	cin >> beta;
	cout << "Inserisci ∆t (s) ";
	cin >> dt;
	
	N_it_max = lround(N_cycl_max*(2.0*(N+1))/dt);
	
	
	double* x = new double[N+2];
	double* v = new double[N+2];
	double* x_new = new double[N+2];
	double* v_new = new double[N+2];
	
	double* C = new double[N_modi+1];
	double* dC = new double[N_modi+1];
	
	double E_n;
	
	
	
	ofstream out("out.txt");
	
	// inizializzazione primo step
	
	// punti fissi (j=0 e j=N+1)
	x[0] = 0;
	x[N+1] = 0;
	v[0] = 0;
	v[N+1] = 0;
	x_new[0] = 0;
	x_new[N+1] = 0;
	v_new[0] = 0;
	v_new[N+1] = 0;
	
	C[0] = 0;
	dC[0] = 0;
	
	
	// coordinate punti mobili
	for(int j=1; j<=N; j++)
	{
		x_new[j] = ampl * sin((M_PI*j)/(N+1));
		v_new[j] = 0;
	}
	
	// coefficienti e calcolo dell'energia
	out << N_it*dt/(2.0*(N+1)) << " "; 	// numero corrispondente di cicli alla frequenza fondamentale
	
	for(int n=1; n<=N_modi; n++)
	{
		C[n] = 0;
		dC[n] = 0;
		E_n = 0;
		
		if(n==1){
			
			for(int j=1; j<=N; j++)
			{
				C[n] += sin((n*M_PI*j)/(N+1))*x_new[j];
				dC[n] += sin((n*M_PI*j)/(N+1))*v_new[j];	
			}
		
		
			C[n] *= sqrt(2.0/(N+1));
			dC[n] *= sqrt(2.0/(N+1));
			
			E_n = 0.5*m*dC[n]*dC[n] + 2*k*C[n]*C[n]*sin(n*M_PI/(2.0*(N+1)));
		}
		
		// stampa valore energia
		out << E_n << " ";
		
	}
	out << endl;
	// qui finisce l'inizializzazione del primo step
	
	
	// ciclo periodico per calcolare tutte le energie
	while(N_it < N_it_max)
	{
		// traslo la griglia avanti di dt
		for(int j=1; j<=N; j++)
		{
			x[j] = x_new[j];
			v[j] = v_new[j];
		}
	
	
		// Velocity Verlet: posizioni 
		for(int j=1; j<=N; j++)
		{
			
			x_new[j] = x[j] + v[j]*dt + 0.5*dt*dt*(
				(k/m)*(x[j-1] - 2*x[j] + x[j+1]) + 
				(alpha/m)*(pow(x[j+1]-x[j], 2.0)-pow(x[j]-x[j-1], 2.0))+ 
				(beta/m)*(pow(x[j+1]-x[j], 3.0)-pow(x[j]-x[j-1], 3.0)));
			
		}
		
		// Velocity Verlet: velocità
		for(int j=1; j<=N; j++)
		{
			
			v_new[j] = v[j] + 0.5*dt*(
				(k/m)*(x[j-1] - 2*x[j] + x[j+1]) + 
				(alpha/m)*(pow(x[j+1]-x[j], 2.0)-pow(x[j]-x[j-1], 2.0))+ 
				(beta/m)*(pow(x[j+1]-x[j], 3.0)-pow(x[j]-x[j-1], 3.0))
				) + 0.5*dt*(
				(k/m)*(x_new[j-1] - 2*x_new[j] + x_new[j+1]) + 
				(alpha/m)*(pow(x_new[j+1]-x_new[j], 2.0)-pow(x_new[j]-x_new[j-1], 2.0))+ 
				(beta/m)*(pow(x_new[j+1]-x_new[j], 3.0)-pow(x_new[j]-x_new[j-1], 3.0)));
		}
		
		
		
		// coefficienti e calcolo dell'energia
		N_it++;
		out << N_it*dt/(2.0*(N+1)) << " "; 	// numero corrispondente di cicli alla frequenza fondamentale
	
		for(int n=1; n<=N_modi; n++)
		{
			C[n] = 0;
			dC[n] = 0;
		
			for(int j=1; j<=N; j++)
			{
				C[n] += sin((n*M_PI*j)/(N+1))*x_new[j];
				dC[n] += sin((n*M_PI*j)/(N+1))*v_new[j];
			}
		
		
			C[n] *= sqrt(2.0/(N+1));
			dC[n] *= sqrt(2.0/(N+1));
		
			E_n = 0.5*m*dC[n]*dC[n] + 2*k*C[n]*C[n]*sin(n*M_PI/(2.0*(N+1)))*sin(n*M_PI/(2.0*(N+1)));
			out << E_n << " ";
		
		}
		out << endl;
		
		
	}
	
	
	delete[] x;
	delete[] v;
	delete[] x_new;
	delete[] v_new;
	delete[] C;
	delete[] dC;
	
	
	
	out.close();
	
	
	return 0;
	
	
}
