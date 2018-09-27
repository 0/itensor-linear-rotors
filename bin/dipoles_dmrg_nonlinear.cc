#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "linrot/operators.h"

#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc != 9) {
        printfln("usage: %s <geom> <R> <N> <l_max> <sweep_table> <dH2_goal> <sweeps_min> <sweeps_max>", argv[0]);

        return 1;
    }

    char *geom_path = argv[1];
    Real R = atof(argv[2]);
    int N = atoi(argv[3]);
    int l_max = atoi(argv[4]);
    auto sweep_table = InputGroup(argv[5], "sweeps");
    Real dH2_goal = atof(argv[6]);
    int sweeps_min = atoi(argv[7]);
    int sweeps_max = atoi(argv[8]);

    double *geom = new double[3*N];

    {
        std::ifstream geom_ifs(geom_path);
        for (auto i : range(N)) {
            geom_ifs >> geom[3*i+0] >> geom[3*i+1] >> geom[3*i+2];
        }
    }

    auto sites = LinearRigidRotor(N, {"l_max", l_max, "m_sym", false});

    auto ampo = AutoMPO(sites);
    // Rotational energy.
    for (auto i : range1(N)) {
        ampo += "L2", i;
    }
    // Potential energy.
    for (auto i : range1(N)) {
        for (auto j : range1(i+1, N)) {
            Real rx = geom[3*(j-1)+0] - geom[3*(i-1)+0];
            Real ry = geom[3*(j-1)+1] - geom[3*(i-1)+1];
            Real rz = geom[3*(j-1)+2] - geom[3*(i-1)+2];
            Real r = std::sqrt(rx*rx + ry*ry + rz*rz);
            Real d = R*r;
            Real k = 1.0/(d*d*d);

            add_operator(ampo, LinearRigidRotorSite::compound_op2("D-D", {"rx", rx/r, "ry", ry/r, "rz", rz/r}), i, j, k);
        }
    }
    auto H = toMPO<IQTensor>(ampo, {"Cutoff=", 1e-40});
    printfln("mean MPO bond dimension: %f", averageM(H));

    delete[] geom;

    auto psi = run_dmrg(sites, N, sweep_table, sweeps_min, sweeps_max, H, dH2_goal);
    run_analysis(psi);

    return 0;
}
