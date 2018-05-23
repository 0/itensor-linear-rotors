#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "linrot/operators.h"

#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc != 7) {
        printfln("usage: %s <R> <N> <field_strength> <l_max> <sweep_table> <N_sweeps>", argv[0]);

        return 1;
    }

    Real R = atof(argv[1]);
    int N = atoi(argv[2]);
    Real field_strength = atof(argv[3]);
    int l_max = atoi(argv[4]);
    auto sweep_table = InputGroup(argv[5], "sweeps");
    int N_sweeps = atoi(argv[6]);

    auto sites = LinearRigidRotor(N, {"l_max", l_max, "lp_sym", false});

    auto ampo = AutoMPO(sites);
    // Rotational energy.
    for (auto i : range1(N)) {
        ampo += "L2", i;
    }
    // Potential energy.
    for (auto i : range1(N)) {
        for (auto j : range1(i+1, N)) {
            Real d = R*(j-i);
#ifdef MINIMAGE
            // Minimum image convention for periodic boundary conditions.
            if (j-i > N/2) {
                d = R*N - d;
            }
#endif // MINIMAGE
            Real k = 1.0/(d*d*d);

            add_operator(ampo, LinearRigidRotorSite::compound_op2("D-D lin"), i, j, k);
        }
    }
#ifdef LINFIELD
    // Linear transverse field.
    for (auto i : range1(N)) {
        add_operator(ampo, LinearRigidRotorSite::compound_op1("x"), i, field_strength*(N-i)/(N-1));
    }
#else // LINFIELD
    // Uniform transverse field.
    for (auto i : range1(N)) {
        add_operator(ampo, LinearRigidRotorSite::compound_op1("x"), i, field_strength);
    }
#endif // LINFIELD
    auto H = toMPO<IQTensor>(ampo, {"Cutoff=", 1e-40});
    printfln("mean MPO bond dimension: %f", averageM(H));

    run_dmrg(sites, ampo, N, N_sweeps, sweep_table, H);

    return 0;
}
