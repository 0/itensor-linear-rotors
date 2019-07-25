#ifndef __LINROT_DMRG_H
#define __LINROT_DMRG_H

#include "itensor/all.h"
#include "linrot/linrigrot.h"

using namespace itensor;


void dmrg_sweep(MPS& psi, MPO const& H, InputGroup& sweep_table, int num_sweeps, int skip_sweeps, std::vector<MPS> ortho_wfs);

void run_analysis(LinearRigidRotor const& sites, MPS& psi);

void dump_coefficients(LinearRigidRotor const& sites, MPS const& psi);

void run_sampling(LinearRigidRotor const& sites, MPS& psi, int num_samples);

MPS embiggen(LinearRigidRotor const& sites1, MPS const& mps1, LinearRigidRotor const& sites2);

#endif // __LINROT_DMRG_H
