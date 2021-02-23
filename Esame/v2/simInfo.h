#ifndef simInfo_h
#define simInfo_h

#include "Singleton.h"

class simInfo: public Singleton<simInfo> {

  friend class Singleton<simInfo>;

 public:
  
  // ask info about simulation via iostream
  void askInfo();
  
  
  // iteration parameters
  double dt() const 		{ return dt_;			}
  int nTimeSteps() const	{ return nTimeSteps_; 	}
  int freqPrint() const		{ return freqPrint_;	}
  
  // simulation parameters
  double qx() const 		{ return qx_;			}
  double qy() const 		{ return qy_;			}
  double V0() const 		{ return V0_;			}
  
  


 private:

  // private constructor and destructor for singleton
  simInfo();
  // deleted copy constructor and assignment to prevent unadvertent copy
  simInfo           ( const simInfo& x ) = delete;
  simInfo& operator=( const simInfo& x ) = delete;

  // destructor
  ~simInfo() override;



  // iteration parameters
  double dt_;				// time interval
  int nTimeSteps_;			// simulation time steps
  int freqPrint_;			// printing frequence
  
  // simulation parameters
  double qx_;
  double qy_;
  double V0_;				// barrier potential
  

};

#endif

