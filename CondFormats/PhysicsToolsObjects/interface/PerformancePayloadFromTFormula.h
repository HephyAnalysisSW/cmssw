#ifndef PerformancePayloadFromTFormula_h
#define PerformancePayloadFromTFormula_h

#include "CondFormats/Serialization/interface/Serializable.h"

#include "CondFormats/PhysicsToolsObjects/interface/PhysicsTFormulaPayload.h"
#include "CondFormats/PhysicsToolsObjects/interface/PerformancePayload.h"

#include <algorithm>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include "TFormula.h"
#include "CondFormats/PhysicsToolsObjects/interface/BinningPointByMap.h"

class PerformancePayloadFromTFormula : public PerformancePayload {
//  class PerformancePayloadFromTFormula : public PerformancePayload, public PhysicsPerformancePayload {
 public:

  static const int InvalidPos;

  PerformancePayloadFromTFormula(const std::vector<PerformanceResult::ResultType>& r, 
				 const std::vector<BinningVariables::BinningVariablesType>& b  ,  
				 PhysicsTFormulaPayload& in) : pl(in), results_(r), variables_(b) {
    initialize();
  }

  PerformancePayloadFromTFormula(){}

  void initialize() override; 

 PerformancePayloadFromTFormula(const PerformancePayloadFromTFormula& b) :
    pl(b.pl), 
    results_(b.results_), 
    variables_(b.variables_), 
    compiledFormulas_(b.compiledFormulas_) {}

  ~PerformancePayloadFromTFormula() override{    
    compiledFormulas_.clear();
  }



  float getResult(PerformanceResult::ResultType,const BinningPointByMap&) const override ; // gets from the full payload
  
  virtual bool isParametrizedInVariable(const BinningVariables::BinningVariablesType p)  const {
    return (limitPos(p) != PerformancePayloadFromTFormula::InvalidPos);
  }
  
  bool isInPayload(PerformanceResult::ResultType,const BinningPointByMap&) const override ;
  
  const PhysicsTFormulaPayload & formulaPayload() const {return pl;}
  
  void printFormula(PerformanceResult::ResultType res) const;
  

 protected:
  
  virtual std::vector<BinningVariables::BinningVariablesType> myBinning() const {return variables_;}

  virtual int limitPos(const BinningVariables::BinningVariablesType b) const {
    std::vector<BinningVariables::BinningVariablesType>::const_iterator p;
    p = find(variables_.begin(), variables_.end(), b);
    if (p == variables_.end()) return PerformancePayloadFromTFormula::InvalidPos;
    return ((p-variables_.begin()));
    
  }

  virtual int resultPos(PerformanceResult::ResultType r) const {
    std::vector<PerformanceResult::ResultType>::const_iterator p;
    p = find (results_.begin(), results_.end(), r);
    if ( p == results_.end()) return PerformancePayloadFromTFormula::InvalidPos;
      return ((p-results_.begin()));
  }


  bool isOk(const BinningPointByMap& p) const; 

  PhysicsTFormulaPayload pl;
  //
  // the variable mapping
  //
  std::vector<PerformanceResult::ResultType> results_;
  std::vector<BinningVariables::BinningVariablesType> variables_;
  //
  // the transient part
  //
  std::vector< std::shared_ptr<const TFormula> > compiledFormulas_ COND_TRANSIENT;

 COND_SERIALIZABLE;
};

#endif
