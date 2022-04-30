#ifndef CMS1703_H_
#define CMS1703_H_
// AUTHOR: Avinash Verma
//  EMAIL: avinash.verma@students.iiserpune.ac.in
#include "AnalysisBase.h"

class Cms1703 : public AnalysisBase {
  public:
    Cms1703() : AnalysisBase()  {}               
    ~Cms1703() {}
  
    void initialize();
    void analyze();        
    void finalize();

  private:
};

#endif
