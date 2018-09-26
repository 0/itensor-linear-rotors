#ifndef __LINROT_DMRG_H
#define __LINROT_DMRG_H

#include "itensor/all.h"
#include "linrot/linrigrot.h"

using namespace itensor;


IQMPS run_dmrg(SiteSet const& sites, int N, InputGroup& sweep_table, int sweeps_min, int sweeps_max, IQMPO const& H, Real dH2_goal);

#endif // __LINROT_DMRG_H
