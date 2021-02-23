/*
	Paolo Zinesi

	Compute time evolution for a 2D-system in a state
	described by a wave function in a box with infinite
	potential at borders (wave function is null at the 
	borders)
	
	Eigen library is needed to run this program
	(https://gitlab.com/libeigen/eigen/-/archive/3.3.9/eigen-3.3.9.zip)
	

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
#include <time.h>

#include "simInfo.h"
#include "matrixUtil.h"

using SpMat = matrixUtil::SpMat;
using T = matrixUtil::T;
using VectComplex = matrixUtil::VectComplex;


// space dimensions
int Nx = 100 + 1;
int Ny = 100 + 1;
double Lx = 500.0;
double Ly = 400.0;
double hx = Lx / (Nx-1);	// ci sono N punti e 
double hy = Ly / (Ny-1);	// N-1 spazi tra di essi

double ax = 250.0;		//
double bx = 260.0;		// dimensioni
double ay = 0.0;			// barriera 2D
double by = Ly;			//

// initial parameters
double x0i = 200.0;
double y0i = 200.0;
double sigmax = 20.0;
double sigmay = 20.0;


int main() {
	
	// iter ≤ nTimeSteps
	int iter = 0;
	
	// create and fill singleton
	// containing infos
	static simInfo* sI = simInfo::instance();
	sI -> askInfo();
	
	// recover variables
	int nTimeSteps = sI -> nTimeSteps();
	int freqPrint  = sI -> freqPrint();
	
	
	// start timer
	std::cout << "Attendere ~ " << (int)(0.087 * nTimeSteps + 4.3) << " secondi" << std::endl;
	time_t timei = time(NULL);
	
	
	// output file
	std::fstream fout;
	fout.open("TDSchr2D.dat",std::ios::out);
	fout.precision(5);
	
	
	// create vector to store results
	VectComplex psi_n(Nx*Ny);
	
	
	// initialize psi_0 with bi-dimensional gaussian 
	// using global variables
	matrixUtil::init2DGaussian(psi_n);
	
	
	// normalize gaussian using Simpson 2D
	matrixUtil::normalize(psi_n);
	
	
	// print on file 
	matrixUtil::filePrint(psi_n, fout);
	
	
	// variables for equation
	// M * psi_n = F 
	SpMat M(Nx*Ny,Nx*Ny);
	VectComplex F(Nx*Ny);
	std::vector<T> tripletList;
	
	
	// fill triplet list
	matrixUtil::fillTripletList(tripletList);
	
	
	// build Spare Matrix
	M.setFromTriplets(tripletList.begin(), tripletList.end());
	
	// build matrix solver
	Eigen::SparseLU<SpMat> solver;
	solver.compute(M);


	
	// begin iterations
	while(++iter < nTimeSteps)
	{
	
		// fill vector F
		matrixUtil::fillF(F, psi_n);
		
	
		// solve the M * psi_n = F equation
		psi_n = solver.solve(F);
		
		
		// print on file
		if(iter % freqPrint == 0) matrixUtil::filePrint(psi_n, fout);
		
	}
	
	
	// stop timer
	time_t timef = time(NULL);
	std::cout << timef-timei << " secondi di esecuzione " << std::endl;
	
	
	fout.close();
	return 0;
	
}