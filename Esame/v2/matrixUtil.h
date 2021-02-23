#ifndef matrixUtil_h
#define matrixUtil_h

#include <Eigen/Sparse>
#include <complex>
#include <fstream>
#include <vector>


class matrixUtil {
	
  public:
		
	  matrixUtil();
	  ~matrixUtil();
	  
	  
	  // typedef for matrix
	  typedef Eigen::SparseMatrix< std::complex<double> >					SpMat;
	  typedef Eigen::Triplet< std::complex<double> > 						T;
	  typedef Eigen::Matrix< std::complex<double>, Eigen::Dynamic, 1 >		VectComplex;
	  
	  
	  // initialize psi with bi-dimensional gaussian
 	  static void init2DGaussian( VectComplex& psi );
	  
	  // normalize 2D graph
	  static void normalize( VectComplex& psi );
	  
	  // print on file
	  static void filePrint( VectComplex& psi, std::fstream& of );
	  
	  // fill triplet list
	  static void fillTripletList( std::vector<T>& tripletList );
	  
	  // fill known vector F
	  // of M * psi = F
	  static void fillF( VectComplex& F , VectComplex& psi );
	  
};

#endif
