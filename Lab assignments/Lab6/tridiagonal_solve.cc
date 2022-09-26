#include <complex>
#include <cmath>

using namespace std;


// funzione che risolve un'equazione di matrice tridiagonale
// -- convenzioni --: lower[0]=0 && upper[n-1]=0

void tridiagonal_solve(int n, complex<double>* upper, complex<double>* lower, 
						complex<double>* diag, complex<double>* sol, complex<double>* b) {
							
							// vettori per contenere i coefficienti
							complex<double>* alpha = new complex<double>[n-1];
							complex<double>* beta = new complex<double>[n-1];
							
							
							// inizializzo il primo termine di alpha e beta
							alpha[0] = (-diag[0]/upper[0]);
							beta[0] = (b[0]/upper[0]);
							
							
							// calcolo restanti termini di alpha e beta
							for(int j=1; j<=n-2; j++)
							{
								alpha[j] = (- diag[j]/upper[j] - lower[j]/(upper[j]*alpha[j-1]));
								beta[j] = (b[j]/upper[j] + (lower[j]*beta[j-1])/(upper[j]*alpha[j-1]));
							}
							
							
							// calcolo primo termine della soluzione
							sol[n-1] = (b[n-1]/lower[n-1] + beta[n-2]/alpha[n-2])/(1.0/alpha[n-2] + diag[n-1]/lower[n-1]);
							
							
							// calcolo resto della soluzione
							for(int j=n-2; j>=0; j--)
							{
								sol[j] = (sol[j+1]-beta[j])/alpha[j];
							}
							
							
							// delete arrays
							delete[] alpha;
							delete[] beta;
							
							
							return;
							
						}							




