/*
Si implementi un codice basato sull’algoritmo di Metropolis per studiare 
la magnetizzazione nel modello di Ising al variare della temperatura.
Si considerino condizioni al contorno periodiche. Ad ogni passo cambiare 
l’orientazione di un solo spin. Considerare N x = N y = 10, 20, 50, J = 1 e T in [0.3].
*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string>

constexpr int N = 50;
constexpr int N_it = 100e+06;



using namespace std;



inline int i_of_l(int l) { return (l-(l%N))/N; }	// i righe della matrice, 0≤i≤N-1
inline int j_of_l(int l) { return (l%N); }			// j colonne della matrice, 0≤j≤N-1
inline int l_of_ij(int i, int j) { return (i*N+j); }

inline void LCG_int (unsigned int &x, unsigned int A, unsigned int M)
{
	x = (A*x) % M;
	return;
}

inline double LCG_01 (unsigned int &x, unsigned int A, unsigned int M)
{
	x = (A*x) % M;
	return (1.0*x)/M;
}



int main()
{
	// parametri LCG
	unsigned int a = pow(7,5);
	unsigned int m = pow(2,31) - 1;
	unsigned int x_i = time(NULL) % m;
	
	
	// informazioni sulla temperatura
	double T;
	cout << "inserisci T ";
	cin >> T;
	double beta = 1.0/(T); /* k_B = 1 */
	
	
	// file di output
	string outfileE, outfileM, outfile_Plot;
	outfileE = "E" + to_string(lround(T*1000)) + "m.dat";
	outfileM = "M" + to_string(lround(T*1000)) + "m.dat";
	outfile_Plot = "sigma" + to_string(lround(T*1000)) + "m.dat";
	
	
	ofstream oE(outfileE);
	ofstream oM(outfileM);
	ofstream oPlot(outfile_Plot);
	oE.precision(6);
	oM.precision(6);
	
	
	int* sigma = new int[N*N];
	int sigma_tmp;
	int N_trial;
	
	
	// inizializza vettore
	int l;
	for(l=0; l<N*N; l++)
	{
		LCG_int(x_i,a,m);				// genero x_i casualmente
		sigma[l] = pow(-1, x_i%2);		// sigma[l] = ±1
		
	}
	
	
	// inizializzo calcolo E,M
	double E_sigma = 0.0, M_sigma = 0.0;
	for(l=0; l<N*N; l++)
	{
		
		M_sigma += sigma[l];
		
		
		// implementazione Periodic-Boundary-Conditions
		
		sigma_tmp = sigma[l_of_ij( (N+i_of_l(l)-1)%N , j_of_l(l) )];		// interazione con l-N
		E_sigma += (-0.5)*sigma[l]*sigma_tmp;
		
		sigma_tmp = sigma[l_of_ij( (N+i_of_l(l)+1)%N , j_of_l(l) )];		// interazione con l+N	
		E_sigma += (-0.5)*sigma[l]*sigma_tmp;
		
		sigma_tmp = sigma[l_of_ij( i_of_l(l) , (N+j_of_l(l)-1)%N )];		// interazione con l-1
		E_sigma += (-0.5)*sigma[l]*sigma_tmp;

		sigma_tmp = sigma[l_of_ij( i_of_l(l) , (N+j_of_l(l)+1)%N )];		// interazione con l+1
		E_sigma += (-0.5)*sigma[l]*sigma_tmp;
		
		
		
	}
	
	double H_sigma = E_sigma;
	E_sigma /= (1.0*N*N);
	M_sigma /= (1.0*N*N);
	
	// fine inizializzazioni
	
	
	// inizio ricorsione
	
	double deltaE, deltaM;
	double deltaH;
	int i=0;
	
	while(i++ <= N_it)
	{
		LCG_int(x_i,a,m);				// genero x_i casualmente
		N_trial = x_i % (N*N);
		
		
		// calcolo delta(H,M,E)
		deltaH = 0.0;
		
		deltaH += 2.0*sigma[N_trial]*sigma[l_of_ij( (N+i_of_l(N_trial)-1)%N , j_of_l(N_trial) )];		// interazione con l-N
		deltaH += 2.0*sigma[N_trial]*sigma[l_of_ij( (N+i_of_l(N_trial)+1)%N , j_of_l(N_trial) )];		// interazione con l+N
		deltaH += 2.0*sigma[N_trial]*sigma[l_of_ij( i_of_l(N_trial) , (N+j_of_l(N_trial)-1)%N )];		// interazione con l-1
		deltaH += 2.0*sigma[N_trial]*sigma[l_of_ij( i_of_l(N_trial) , (N+j_of_l(N_trial)+1)%N )];		// interazione con l+1
		
		
		deltaE = 1.0*deltaH/(N*N);
		
		deltaM = -2.0*sigma[N_trial]/(N*N);
		
		
		// implementazione Metropolis
		
		if (deltaH <= 0)		// accettazione sigma_trial
		{
			sigma[N_trial] *= (-1);
			H_sigma += deltaH;
			E_sigma += deltaE;
			M_sigma += deltaM;
			
		}
		else if (exp(-beta*deltaH) >= LCG_01(x_i,a,m))		// accettazione sigma_trial
		{
			sigma[N_trial] *= (-1);
			H_sigma += deltaH;
			E_sigma += deltaE;
			M_sigma += deltaM;
			
		}
		else		// rifiuto di sigma_trial
		{
			sigma[N_trial] *= (1);
			H_sigma += 0.0;
			E_sigma += 0.0;
			M_sigma += 0.0;
			
		}
		
		
		// stampa
		if(i%1000 == 0)
		{
			oE << i << " " << E_sigma << endl;
			oM << i << " " << M_sigma << endl;
			
		}
		
		if(i%1000 == 0 && i < 1.0e+06)
		{
			// plot matrix of sigmas
			for(int i_l=0; i_l<N; i_l++)
			{
				for(int j_l=0; j_l<N; j_l++)
				{
					oPlot << sigma[l_of_ij(i_l, j_l)] << " ";
				}
				oPlot << endl;
			}
			oPlot << endl << endl;
			
		}
		
		
	}
	
	

	oE.close();
	oM.close();
	oPlot.close();

	delete[] sigma;
	return 0;
	
}