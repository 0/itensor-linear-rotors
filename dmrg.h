#ifndef __LINROT_DMRG_H
#define __LINROT_DMRG_H

#include "itensor/all.h"
#include "linrot/linrigrot.h"

using namespace itensor;


void run_dmrg(SiteSet const& sites, AutoMPO const& ampo, int N, int N_sweeps, InputGroup& sweep_table, IQMPO const& H);

#endif // __LINROT_DMRG_H
