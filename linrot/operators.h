#ifndef __LINROT_OPERATORS_H
#define __LINROT_OPERATORS_H

#include "itensor/mps/autompo.h"


struct OneSiteTerm {
    itensor::Cplx k;
    std::string op;
};

class OneSiteOperator {
    std::vector<OneSiteTerm> terms;

public:
    typedef std::vector<OneSiteTerm>::const_iterator const_iterator;

    OneSiteOperator() { }

    void add_term(itensor::Cplx k, std::string op) {
        terms.push_back({k, op});
    }

    const_iterator begin() const { return terms.begin(); }
    const_iterator end() const { return terms.end(); }
};


struct TwoSiteTerm {
    itensor::Cplx k;
    OneSiteOperator op1, op2;
};

class TwoSiteOperator {
    OneSiteOperator (*op_gen)(std::string const&);
    std::vector<TwoSiteTerm> terms;

public:
    typedef std::vector<TwoSiteTerm>::const_iterator const_iterator;

    TwoSiteOperator(OneSiteOperator (*op_gen)(std::string const&)) : op_gen(op_gen) { }

    void add_term(itensor::Cplx k, std::string op1, std::string op2) {
        terms.push_back({k, op_gen(op1), op_gen(op2)});
    }

    const_iterator begin() const { return terms.begin(); }
    const_iterator end() const { return terms.end(); }
};


void add_operator(itensor::AutoMPO& ampo, OneSiteOperator op, int i, itensor::Cplx k);

void add_operator(itensor::AutoMPO& ampo, TwoSiteOperator op, int i, int j, itensor::Cplx k);

#endif // __LINROT_OPERATORS_H
