#include "dmrg.h"


void run_dmrg(SiteSet const& sites, AutoMPO const& ampo, int N, int N_sweeps, InputGroup& sweep_table, IQMPO const& H) {
    auto state = InitState(sites);
    for (auto i : range1(N)) {
        state.set(i, "l0m0");
    }
    auto psi = IQMPS(state);

    auto sweeps = Sweeps(N_sweeps, sweep_table);
    println();
    println(sweeps);

    auto energy = dmrg(psi, H, sweeps, "Quiet");

    println();
    printfln("E0 = %.15f", energy);

    auto dH2 = N*(overlap(psi, H, H, psi)/(energy*energy) - 1.0);
    printfln("dH2 = %.15f", dH2);

    // Entanglement entropy.
    psi.position(N/2);
    auto wf = psi.A(N/2)*psi.A(N/2+1);
    auto U = psi.A(N/2);

    IQTensor S, V;
    auto spectrum = svd(wf, U, S, V);

    Real SvN = 0.0;
    Real S2 = 0.0;
    for(auto p : spectrum.eigs()) {
        if(p < 1e-12) {
            continue;
        }

        SvN += -p*log(p);
        S2 += p*p;
    }
    S2 = -log(S2);

    printfln("SvN = %.15f", SvN);
    printfln("S2 = %.15f", S2);

    // Orientational correlation.
    auto OC_ampo = AutoMPO(sites);
    for (auto i : range1(N)) {
        for (auto j : range1(i+1, N)) {
            add_operator(OC_ampo, LinearRigidRotorSite::compound_op2("dot product"), i, j, 1.0);
        }
    }
    auto OC = toMPO<IQTensor>(OC_ampo, {"Cutoff=", 1e-40});

    printfln("orientational correlation = %.15f", overlap(psi, OC, psi)*2.0/(N*(N-1)));

    // Spatial correlations.
    Real *corr = new Real[N*N];

    for (auto i : range1(N)) {
        for (auto j : range1(i, N)) {
            corr[N*(i-1)+(j-1)] = 0.0;
        }
    }

    for (const auto& term_i : LinearRigidRotorSite::compound_op1("z")) {
        for (const auto& term_j : LinearRigidRotorSite::compound_op1("z")) {
            for (auto i : range1(N)) {
                psi.position(i);

                auto op_i = sites.op(term_i.op, i);

                // Diagonal.
                {
                    auto op_i1 = op_i;
                    auto op_i2 = sites.op(term_j.op, i);

                    auto C = psi.A(i) * op_i1 * prime(op_i2) * dag(prime(psi.A(i), Site, 2));
                    auto zz = term_i.k * term_j.k * C.cplx();

                    if (std::abs(zz.imag()) > 1e-12) {
                        println("WARNING: Complex correlation!");
                    }

                    corr[(N+1)*(i-1)] += zz.real();
                }

                // Off-diagonal.
                auto idx = commonIndex(psi.A(i), psi.A(i+1), Link);
                auto C = psi.A(i) * op_i * dag(prime(psi.A(i), Site, idx));

                for (auto j : range1(i+1, N)) {
                    auto op_j = sites.op(term_j.op, j);

                    if (j-1 > i) {
                        C *= psi.A(j-1);
                        C *= dag(prime(psi.A(j-1), Link));
                    }

                    auto D = C * psi.A(j);
                    D *= op_j;
                    auto idx = commonIndex(psi.A(j), psi.A(j-1), Link);
                    D *= dag(prime(psi.A(j), idx, Site));
                    auto zz = term_i.k * term_j.k * D.cplx();

                    if (std::abs(zz.imag()) > 1e-12) {
                        println("WARNING: Complex correlation!");
                    }

                    corr[N*(i-1)+(j-1)] += zz.real();
                }
            }
        }
    }

    for (auto i : range1(N)) {
        for (auto j : range1(i, N)) {
            printfln("<z%04d z%04d> = %.15f", i, j, corr[N*(i-1)+(j-1)]);
        }
    }

    delete[] corr;
}
