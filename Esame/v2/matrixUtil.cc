#include "matrixUtil.h"
#include <math.h>
#include "simInfo.h"


// space dimensions
extern int Nx;
extern int Ny;
extern double hx;			// ci sono N punti e 
extern double hy;			// N-1 spazi tra di essi

extern double ax;			//
extern double bx;			// dimensioni
extern double ay;			// barriera 2D
extern double by;			//

// initial parameters
extern double x0i;
extern double y0i	;
extern double sigmax;
extern double sigmay;


// inline definitions (for multi-index)
inline int i_of_l (int l) 		 	{ 	return l % Nx; 				}		// 0 ≤ i ≤ Nx-1
inline int j_of_l (int l) 		 	{ 	return (l - (l % Nx)) / Nx; }		// 0 ≤ j ≤ Ny-1
inline int l_of_ij (int i, int j) 	{ 	return Nx * j + i; 			}		// 0 ≤ l ≤ Nx*Ny-1


matrixUtil::matrixUtil(){
}

matrixUtil::~matrixUtil() {
}

// initialize psi with bi-dimensional gaussian 
void matrixUtil::init2DGaussian( VectComplex& psi ) {
	
	static simInfo* sI = simInfo::instance();
	double qx = sI -> qx();
	double qy = sI -> qy();	
	
	
	int l;
	for(l=0; l<Nx*Ny; l++)
	{	
		// intern
		if(i_of_l(l) != 0 && i_of_l(l) != Nx-1 && j_of_l(l) != 0 && j_of_l(l) != Ny-1)
			psi[l] =	std::polar(1.0, qx * i_of_l(l) * hx + qy * j_of_l(l) * hy) *
						exp(-0.5 * pow((i_of_l(l) * hx - x0i) / sigmax, 2.0)) * 
						exp(-0.5 * pow((j_of_l(l) * hy - y0i) / sigmay, 2.0));
		
		// border
		else psi[l] = std::complex<double>(0.0, 0.0);
		
	}
	
	return;
	
}


// normalize 2D graph
void matrixUtil::normalize( VectComplex& psi ) {
	
	int l;
	double norm = 0.0;
	
	// compute norm
	for(l=0; l<Nx*Ny; l++)
	{	
		if 		(i_of_l(l) % 2 == 0 && j_of_l(l) % 2 == 0) 	norm += 	std::norm(psi[l]);
		else if (i_of_l(l) % 2 == 1 && j_of_l(l) % 2 == 1) 	norm += 4 * std::norm(psi[l]);
		else 												norm += 2 * std::norm(psi[l]);
		
	}
	norm *= (4.0 * hx * hy / 9.0);
	norm = sqrt(norm);
	
	
	// normalize
	for(l=0; l<Nx*Ny; l++) psi[l] /= norm;
	
	
	return;
	
}


// print on file
void matrixUtil::filePrint( VectComplex& psi, std::fstream& of ) {
	
	int i,j;
	
	for(j=0; j<Ny; j++)
	{
		for(i=0; i<Nx; i++)
		{
			of << std::abs(psi[l_of_ij(i,j)]) << " ";
		}
		of << std::endl;
	}
	of << std::endl << std::endl;
	
	
	return;
	
}


// fill triplet list
void matrixUtil::fillTripletList( std::vector<T>& tripletList ) {
	
	tripletList.clear();
	tripletList.reserve(5*Nx*Ny);
	
	int l;
	std::complex<double> dij;
	
	static simInfo* sI = simInfo::instance();
	double dt = sI -> dt();
	double V0 = sI -> V0();
	
	
	// fill list according to potential
	for(l=0; l<Nx*Ny; l++)
	{
		// intern
		if(i_of_l(l) != 0 && i_of_l(l) != Nx-1 && j_of_l(l) != 0 && j_of_l(l) != Ny-1)
		{
			
			dij = std::complex<double>(0.0, -4.0 / dt) + 2.0/(hx*hx) + 2.0/(hy*hy) ;
			
			// Vij ≠ 0 only for ax ≤ x ≤ bx && ay ≤ y ≤ by
			if(	ax <= i_of_l(l) * hx && i_of_l(l) * hx <= bx &&
				ay <= j_of_l(l) * hy && j_of_l(l) * hy <= by 	) dij += 2.0 * V0;
			
			
			// diagonal
			tripletList.push_back(T(l, l, dij ));
			
			// non-diagonal
			tripletList.push_back(T(l, l_of_ij( i_of_l(l) + 1, j_of_l(l) ), -1.0/(hx*hx) ));
			tripletList.push_back(T(l, l_of_ij( i_of_l(l) - 1, j_of_l(l) ), -1.0/(hx*hx) ));
			tripletList.push_back(T(l, l_of_ij( i_of_l(l), j_of_l(l) + 1 ), -1.0/(hy*hy) ));
			tripletList.push_back(T(l, l_of_ij( i_of_l(l), j_of_l(l) - 1 ), -1.0/(hy*hy) ));
			
			
		}	
		else // border
			tripletList.push_back(T(l, l, 1.0 ));
	
	}
	
	return;
	
}


// fill known vector F
// of M * psi = F
void matrixUtil::fillF( VectComplex& F , VectComplex& psi ) {
	
	int l;
	
	static simInfo* sI = simInfo::instance();
	double dt = sI -> dt();
	double V0 = sI -> V0();
	
	
	// fill F
	for(l=0; l<Nx*Ny; l++)
	{
		// intern
		if(i_of_l(l) != 0 && i_of_l(l) != Nx-1 && j_of_l(l) != 0 && j_of_l(l) != Ny-1)
		{
			
			// Vij ≠ 0 only for ax ≤ x ≤ bx && ay ≤ y ≤ by
			F[l] = std::complex<double>(0.0, -4.0 / dt ) - 2.0/(hx*hx) - 2.0/(hy*hy);
			
			if(	ax <= i_of_l(l) * hx && i_of_l(l) * hx <= bx &&
				ay <= j_of_l(l) * hy && j_of_l(l) * hy <= by 	) F[l] += (-2.0 * V0);
			
			// diagonal
			F[l] *= psi[l];
			
			// non-diagonal
			F[l] += psi[l_of_ij( i_of_l(l) + 1 , j_of_l(l) )] * 1.0 / (hx*hx);
			F[l] += psi[l_of_ij( i_of_l(l) - 1 , j_of_l(l) )] * 1.0 / (hx*hx);
			F[l] += psi[l_of_ij( i_of_l(l) , j_of_l(l) + 1 )] * 1.0 / (hy*hy);
			F[l] += psi[l_of_ij( i_of_l(l) , j_of_l(l) - 1 )] * 1.0 / (hy*hy);
			
							
		}
		else // border
			F[l] = 0.0;
			
	}
	
	return;
	
}