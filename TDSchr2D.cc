/*
	Paolo Zinesi

	Compute time evolution for a 2D-system in a state
	described by a wave function in a box with infinite
	potential at borders (wave function is null at the 
	borders)


	Conventions:
	
	(atomic unities)
	h/(2*M_PI) = 1
	m = 1
	
*/
	

#include <iostream>
#include <vector>
#include <Eigen/Sparse>
#include <complex>
#include <math.h>
#include <fstream>
#include <time.h>

// space dimensions
constexpr int Nx = 100 + 1;
constexpr int Ny = 100 + 1;
constexpr double Lx = 500.0;
constexpr double Ly = 400.0;
constexpr double hx = Lx / (Nx-1);	// ci sono N punti e 
constexpr double hy = Ly / (Ny-1);	// N-1 spazi tra di essi
constexpr double a = 250.0;
constexpr double b = 260.0;

// initial parameters
constexpr double x0i = 200.0;
constexpr double y0i = 200.0;
constexpr double sigmax = 20.0;
constexpr double sigmay = 20.0;

// inline definitions (for multi-index)
inline int i_of_l (int l) 		 	{ 	return l % Nx; 				}		// 0 ≤ i ≤ Nx-1
inline int j_of_l (int l) 		 	{ 	return (l - (l % Nx)) / Nx; }		// 0 ≤ j ≤ Ny-1
inline int l_of_ij (int i, int j) 	{ 	return Nx * j + i; 			}		// 0 ≤ l ≤ Nx*Ny-1


// typedef for matrix
typedef Eigen::SparseMatrix< std::complex<double> >					SpMat;
typedef Eigen::Triplet< std::complex<double> > 						T;
typedef Eigen::Matrix< std::complex<double>, Eigen::Dynamic, 1 >	VectComplex;

int main() {
	
	// iter ≤ nTimeSteps
	int iter = 0;
	int i,j,l;
	
	
	// simulation parameters
	double dt;
	int nTimeSteps, freqPrint;
	std::cout << "dt (default 1.0) ";
	std::cin >> dt;
	dt = (dt ? dt : 1.0 );
	std::cout << "nTimeSteps (default 500) ";
	std::cin >> nTimeSteps;
	nTimeSteps = (nTimeSteps ? nTimeSteps : 500 );
	std::cout << "freqPrint (default 5) ";
	std::cin >> freqPrint;
	freqPrint = (freqPrint ? freqPrint : 5 );
	
	double qx, qy;
	std::cout << "qx in [ " << -M_PI/hx << " , " << M_PI/hx << "] ";
	std::cin >> qx;
	std::cout << "qy in [ " << -M_PI/hy << " , " << M_PI/hy << "] ";
	std::cin >> qy;
	
	double V0;
	std::cout << "V0 (soglia ~ " << 0.5/(hx*hy) << ") ";
	std::cin >> V0;
	std::cout << "Attendere ~ " << (int)(0.087 * nTimeSteps + 4.3) << " secondi" << std::endl;
	time_t timei = time(NULL);
	
	
	// output file
	std::fstream fout;
	fout.open("TDSchr2D.dat",std::ios::out);
	fout.precision(5);
	
	
	// create vectors to store results
	VectComplex psi_n(Nx*Ny);
	
	
	// initialize psi_0 with bi-dimensional gaussian 
	for(l=0; l<Nx*Ny; l++)
	{	
		// intern
		if(i_of_l(l) != 0 && i_of_l(l) != Nx-1 && j_of_l(l) != 0 && j_of_l(l) != Ny-1)
			psi_n[l] =	std::polar(1.0, qx * i_of_l(l) * hx + qy * j_of_l(l) * hy) *
						exp(-0.5 * pow((i_of_l(l) * hx - x0i) / sigmax, 2.0)) * 
						exp(-0.5 * pow((j_of_l(l) * hy - y0i) / sigmay, 2.0));
		
		// border
		else psi_n[l] = std::complex<double>(0.0, 0.0);
		
	}
	
	// compute norm (special Simpson 2D) and normalize gaussian
	double norm = 0.0;
	std::cout << "norma teorica " << sqrt(M_PI * sigmax * sigmay) << std::endl;
	
	for(l=0; l<Nx*Ny; l++)
	{	
		if 		(i_of_l(l) % 2 == 0 && j_of_l(l) % 2 == 0) 	norm += 	std::norm(psi_n[l]);
		else if (i_of_l(l) % 2 == 1 && j_of_l(l) % 2 == 1) 	norm += 4 * std::norm(psi_n[l]);
		else 												norm += 2 * std::norm(psi_n[l]);
		
	}
	norm *= (4.0 * hx * hy / 9.0);
	norm = sqrt(norm);
	std::cout << "norma calcolata iniziale " << norm << std::endl;
	for(l=0; l<Nx*Ny; l++) psi_n[l] /= norm;
	
	
	
	
	// print on file 
	for(j=0; j<Ny; j++)
	{
		for(i=0; i<Nx; i++)
		{
			fout << std::abs(psi_n[l_of_ij(i,j)]) << " ";
		}
		fout << std::endl;
	}
	fout << std::endl << std::endl;
	
	
	// utility parameters for equation
	// -------------------------
	// ----- M * psi_n = F -----
	// -------------------------
	
	SpMat M(Nx*Ny,Nx*Ny);
	VectComplex F(Nx*Ny);
	
	std::vector<T> tripletList;
	tripletList.reserve(5*Nx*Ny);
	
	std::complex<double> dij;
	
	
	// fill triplet list
	for(l=0; l<Nx*Ny; l++)
	{
		// intern
		if(i_of_l(l) != 0 && i_of_l(l) != Nx-1 && j_of_l(l) != 0 && j_of_l(l) != Ny-1)
		{
			
			dij = std::complex<double>(0.0, -4.0 / dt) + 2.0/(hx*hx) + 2.0/(hy*hy) ;
			
			// Vij ≠ 0 only for a ≤ x ≤ b
			if(a <= i_of_l(l) * hx && i_of_l(l) * hx <= b ) dij += 2.0 * V0;
			
			
			tripletList.push_back(T(l, l, dij ));
			tripletList.push_back(T(l, l_of_ij( i_of_l(l) + 1, j_of_l(l) ), -1.0/(hx*hx) ));
			tripletList.push_back(T(l, l_of_ij( i_of_l(l) - 1, j_of_l(l) ), -1.0/(hx*hx) ));
			tripletList.push_back(T(l, l_of_ij( i_of_l(l), j_of_l(l) + 1 ), -1.0/(hy*hy) ));
			tripletList.push_back(T(l, l_of_ij( i_of_l(l), j_of_l(l) - 1 ), -1.0/(hy*hy) ));
			
			
		}	
		else // border
			tripletList.push_back(T(l, l, 1.0 ));
	
	}
	
	
	// build Spare Matrix
	M.setFromTriplets(tripletList.begin(), tripletList.end());
	
	// build matrix solver
	Eigen::SparseLU<SpMat> solver;
	solver.compute(M);

	
	// begin iterations
	while(++iter < nTimeSteps)
	{
	
		// create vector F
		for(l=0; l<Nx*Ny; l++)
		{
			// intern
			if(i_of_l(l) != 0 && i_of_l(l) != Nx-1 && j_of_l(l) != 0 && j_of_l(l) != Ny-1)
			{
				
				// Vij ≠ 0 only for a ≤ x ≤ b
				F[l] = std::complex<double>(0.0, -4.0 / dt ) - 2.0/(hx*hx) - 2.0/(hy*hy);
				if(a <= i_of_l(l) * hx && i_of_l(l) * hx <= b ) F[l] += (-2.0 * V0);
				F[l] *= psi_n[l];
				
				
				F[l] += psi_n[l_of_ij( i_of_l(l) + 1 , j_of_l(l) )] * 1.0 / (hx*hx);
				F[l] += psi_n[l_of_ij( i_of_l(l) - 1 , j_of_l(l) )] * 1.0 / (hx*hx);
				F[l] += psi_n[l_of_ij( i_of_l(l) , j_of_l(l) + 1 )] * 1.0 / (hy*hy);
				F[l] += psi_n[l_of_ij( i_of_l(l) , j_of_l(l) - 1 )] * 1.0 / (hy*hy);
								
			}
			else // border
				F[l] = 0.0;
				
		}
	
		// solve the M * psi_n = F equation
		psi_n = solver.solve(F);
		
		
		
		// print on file
		if(iter % 5 != 0) continue;
		for(j=0; j<Ny; j++)
		{
			for(i=0; i<Nx; i++)
			{
				fout << std::abs(psi_n[l_of_ij(i,j)]) << " ";
			}
			fout << std::endl;
		}
		fout << std::endl << std::endl;
		
	}
	// end iterations
	
	
	// compute final norm (special Simpson 2D)
	norm = 0.0;
	for(l=0; l<Nx*Ny; l++)
	{	
		if 		(i_of_l(l) % 2 == 0 && j_of_l(l) % 2 == 0) 	norm += 	std::norm(psi_n[l]);
		else if (i_of_l(l) % 2 == 1 && j_of_l(l) % 2 == 1) 	norm += 4 * std::norm(psi_n[l]);
		else 												norm += 2 * std::norm(psi_n[l]);
		
	}
	norm *= (4.0 * hx * hy / 9.0);
	norm = sqrt(norm);
	std::cout << "norma calcolata fine " << norm << std::endl;
	
	
	// compute total execution time
	time_t timef = time(NULL);
	std::cout << timef-timei << " secondi di esecuzione " << std::endl;
	
	
	fout.close();
	return 0;
	
}