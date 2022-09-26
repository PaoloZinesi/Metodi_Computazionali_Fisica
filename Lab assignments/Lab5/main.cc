#include <iostream>
#include <math.h>
#include <fstream>

// parametri standard per confronto
const double Lx = 10.0, Ly = 10.0;
const int Nx = 100, Ny = 100;
const double px = 4.0, py = 4.0;
const double Lc = 2.0;
const double rho0 = 80.0;


// funzioni per le trasformazioni (i,j)<-->(m)
inline int i_of (int m) { return (m+1-Nx*floor((m*1.0)/Nx)); }		// i(m)
inline int j_of (int m) { return (1+floor((m*1.0)/Nx)); }			// j(m)
inline int m_of_ij(int i, int j) { return (Nx*(j-1)+i-1); }			// m(i,j)


int main()
{
	
	// variabili operative
	double ib = (px/Lx)*Nx;				// x : bordo box di sinistra
	double ie = ((px+Lc)/Lx)*Nx;		// x : bordo box di destra
	double jb = (py/Ly)*Ny;				// y : bordo box inferiore
	double je = ((py+Lc)/Ly)*Ny;		// y : bordo box superiore
	
	// altre variabili operative 
	double hx = Lx/Nx;
	double hy = Ly/Ny;
	double rho_tmp;
	double diff_phi;
	
	
	// distribuzione di rho da considerare
	enum Distr {monop, dip, quadp};
	Distr type;
	int t_buff;
	
	// inserimento distribuzione
	std::cout << "Considerare monopolo (0), dipolo (1) o quadrupolo (2) ";
	std::cin >> t_buff;
	if(t_buff == 0) type = monop;
	else if(t_buff == 1) type = dip;
	else if(t_buff == 2) type = quadp;
	else {
		std::cout << "Non valido ";
		return 0;
	}
	
	// vettore dei risultati
	double* phi = new double[Nx*Ny];
	double* phi_new = new double[Nx*Ny];
	double* phi_0 = new double[Nx*Ny];		// vettore di base
	
	
	
	// inizializzazione di phi (iterazione numero 0)
	int N_it = 0;
	for(int m=0; m<Nx*Ny; m++)
	{
		// calcolo distribuzione
		if(type == Distr::monop)
		{
			// monopolo
			if(ib<=i_of(m) && i_of(m)<=ie && jb<=j_of(m) && j_of(m)<=je) rho_tmp = rho0;
			else rho_tmp = 0;
		}
		else if(type == Distr::dip)
		{
			// dipolo
			if(ib<=i_of(m) && i_of(m)<0.5*(ib+ie) && jb<=j_of(m) && j_of(m)<=je) rho_tmp = -rho0;
			else if(0.5*(ib+ie)<i_of(m) && i_of(m)<=ie && jb<=j_of(m) && j_of(m)<=je) rho_tmp = rho0;
			else rho_tmp = 0;
		}
		else if(type == Distr::quadp)
		{
			// quadrupolo
			if(ib<=i_of(m) && i_of(m)<0.5*(ib+ie) && jb<=j_of(m) && j_of(m)<0.5*(jb+je)) rho_tmp = -rho0;
			else if(0.5*(ib+ie)<i_of(m) && i_of(m)<=ie && 0.5*(jb+je)<j_of(m) && j_of(m)<=je) rho_tmp = -rho0;
			else if(ib<=i_of(m) && i_of(m)<0.5*(ib+ie) && 0.5*(jb+je)<j_of(m) && j_of(m)<=je) rho_tmp = rho0;
			else if(0.5*(ib+ie)<i_of(m) && i_of(m)<=ie && jb<=j_of(m) && j_of(m)<0.5*(jb+je)) rho_tmp = rho0;
			else rho_tmp = 0;
		}
		
		
		phi_0[m] = rho_tmp*((0.5*hx*hx*hy*hy)/(hx*hx+hy*hy));
		phi_new[m] = phi_0[m];
	}
	
	
	// iterazioni successive (finchè non converge)
	do
	{
		
		N_it++;
		
		// libero spazio per la nuova iterazione
		for(int l=0; l<Nx*Ny; l++)
		{
			phi[l] = phi_new[l];
		}
		
		
		// calcolo per ogni componente del vettore
		for(int l=0; l<Nx*Ny; l++)
		{
			phi_new[l] = 0;
			
			
			// se l è sul bordo si evitano calcoli inutili
			if(i_of(l)!=0 && i_of(l)!=Nx-1 && j_of(l)!=0 && j_of(l)!=Ny-1)
			{
				
				// prodotto matriciale M*phi semplificato
				
				//					phi_new[l] += M[l][m] * phi[m]	negli unici casi in cui M[l][m]≠0
				if(i_of(l)+1 <= Nx)	phi_new[l] += ((0.5*hy*hy)/(hx*hx+hy*hy))*phi[m_of_ij(i_of(l)+1, j_of(l))];
				if(1 <= i_of(l)-1)	phi_new[l] += ((0.5*hy*hy)/(hx*hx+hy*hy))*phi[m_of_ij(i_of(l)-1, j_of(l))];
				
				if(j_of(l)+1 <= Ny)	phi_new[l] += ((0.5*hx*hx)/(hx*hx+hy*hy))*phi[m_of_ij(i_of(l), j_of(l)+1)];
				if(1 <= j_of(l)-1)	phi_new[l] += ((0.5*hx*hx)/(hx*hx+hy*hy))*phi[m_of_ij(i_of(l), j_of(l)-1)];				


			
				// somma vettore di inizializzazione
				phi_new[l] += phi_0[l];
					
			}
			
		}
		
		
		// calcolo la norma del vettore differenza
		diff_phi = 0.0;
		for(int l=0; l<Nx*Ny; l++)
		{
			diff_phi += (phi_new[l]-phi[l])*(phi_new[l]-phi[l]);
		}
		
	}
	while (diff_phi > 1e-10);
	
	
	// apri file di output
	std::ofstream out("out.txt");
	
	
	// output risultati
	for(int i=0; i<Nx; i++)
	{
		for(int j=0; j<Ny; j++)
		{
			out << hx*(i+1) << " " << hy*(j+1) << " " << phi_new[m_of_ij(i,j)] << std::endl;
			
		}
		out << std::endl;
	}
	
	// chiudi file di output
	out.close();
	
	
	// elimina array dinamici
	delete[] phi;
	delete[] phi_new;
	delete[] phi_0;
	
	
	return 0;
	
}