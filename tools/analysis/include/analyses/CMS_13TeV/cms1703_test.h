#ifndef CMS1703_TEST_H_
#define CMS1703_TEST_H_
// AUTHOR: Avinash Verma
//  EMAIL: avinash.verma@students.iiserpune.ac.in
#include "AnalysisBase.h"

class Cms1703_test : public AnalysisBase {
  public:
    Cms1703_test() : AnalysisBase()  {}               
    ~Cms1703_test() {}
  
    void initialize();
    void analyze();        
    void finalize();

  private:
    int n=0;
};


class Cms1703_test_CR : public AnalysisBase {
  public:
    Cms1703_test_CR() : AnalysisBase()  {}               
    ~Cms1703_test_CR() {}
    
    void initialize();
    void analyze();        
    void finalize();

  private:
    int n=0;
};

#endif
