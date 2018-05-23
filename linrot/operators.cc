#include "operators.h"


void add_operator(itensor::AutoMPO& ampo, OneSiteOperator op, int i, itensor::Cplx k) {
    for (const auto& term : op) {
        ampo += k*term.k, term.op, i;
    }
}

void add_operator(itensor::AutoMPO& ampo, TwoSiteOperator op, int i, int j, itensor::Cplx k) {
    for (const auto& term : op) {
        for (const auto& term1 : term.op1) {
            for (const auto& term2 : term.op2) {
                ampo += k*term.k*term1.k*term2.k, term1.op, i, term2.op, j;
            }
        }
    }
}
