#ifndef TEST_H_
#define TEST_H_
// AUTHOR: test
//  EMAIL: test
#include "AnalysisBase.h"

class Test : public AnalysisBase {
  public:
    Test() : AnalysisBase()  {}               
    ~Test() {}
  
    void initialize();
    void analyze();        
    void finalize();

  private:
};

#endif
