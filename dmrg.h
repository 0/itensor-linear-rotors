#ifndef __LINROT_DMRG_H
#define __LINROT_DMRG_H

#include "itensor/all.h"
#include "linrot/linrigrot.h"

using namespace itensor;


void dmrg_sweep(IQMPS& psi, IQMPO const& H, InputGroup& sweep_table, int num_sweeps, int skip_sweeps);

void run_analysis(IQMPS& psi);

void dump_coefficients(int l_max, IQMPS const& psi);

void run_sampling(int l_max, IQMPS& psi, int num_samples);

IQMPS embiggen(LinearRigidRotor const& sites1, IQMPS const& mps1, LinearRigidRotor const& sites2);

#endif // __LINROT_DMRG_H
