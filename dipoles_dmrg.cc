#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "linrot/operators.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc != 5) {
        printfln("usage: %s <R> <N> <l_max> <N_sweeps>", argv[0]);

        return 1;
    }

    Real R = atof(argv[1]);
    int N = atoi(argv[2]);
    int l_max = atoi(argv[3]);
    int N_sweeps = atoi(argv[4]);

    auto sites = LinearRigidRotor(N, {"l_max", l_max});

    auto ampo = AutoMPO(sites);
    // Rotational energy.
    for (auto i : range1(N)) {
        ampo += "L2", i;
    }
    // Potential energy.
    for (auto i : range1(N)) {
        for (auto j : range1(i+1, N)) {
            Real d = R*(j-i);
            Real k = 1.0/(d*d*d);

            add_operator(ampo, LinearRigidRotorSite::compound_op2("D-D lin"), i, j, k);
        }
    }
    auto H = toMPO<IQTensor>(ampo, {"Cutoff=", 1e-40});
    printfln("mean MPO bond dimension: %f", averageM(H));

    auto state = InitState(sites);
    for (auto i : range1(N)) {
        state.set(i, "l0m0");
    }
    auto psi = IQMPS(state);

    auto sweeps = Sweeps(N_sweeps);
    sweeps.maxm() = 10, 20, 100, 100, 200, 400, 800;
    sweeps.cutoff() = 1e-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1e-7, 1e-8, 0.0;
    println();
    println(sweeps);

    auto energy = dmrg(psi, H, sweeps, "Quiet");

    println();
    printfln("E0 = %.15f", energy);

    // Orientational correlation.
    auto OC_ampo = AutoMPO(sites);
    for (auto i : range1(N)) {
        for (auto j : range1(i+1, N)) {
            add_operator(OC_ampo, LinearRigidRotorSite::compound_op2("dot product"), i, j, 1.0);
        }
    }
    auto OC = toMPO<IQTensor>(OC_ampo, {"Cutoff=", 1e-40});

    printfln("orientational correlation = %.15f", overlap(psi, OC, psi)*2.0/(N*(N-1)));

    return 0;
}
