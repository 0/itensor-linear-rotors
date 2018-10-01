#ifndef __LINROT_DMRG_H
#define __LINROT_DMRG_H

#include "itensor/all.h"
#include "linrot/linrigrot.h"

using namespace itensor;


IQMPS run_dmrg(SiteSet const& sites, int N, InputGroup& sweep_table, int sweeps_min, int sweeps_max, IQMPO const& H, Real dH2_goal);

void dump_probabilities(IQMPS const& psi, int l_max);

void run_sampling(IQMPS& psi, int l_max, int num_samples);

#endif // __LINROT_DMRG_H
