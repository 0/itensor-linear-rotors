#include "observer.h"
#include "dmrg.h"


Sweeps make_sweeps(InputGroup& sweep_table, int num_sweeps, int skip_sweeps) {
    if (!sweep_table.GotoGroup()) {
        Error("Table \"" + sweep_table.name() + "\" not found");
    }

    // Skip the header.
    sweep_table.SkipLine();

    auto sweeps = Sweeps(num_sweeps);

    int sw_phys = 0;
    int maxm_last, minm_last, niter_last;
    Real cutoff_last, noise_last;

    for (int sw = 1; sw <= num_sweeps; sw++) {
        int maxm_cur, minm_cur, niter_cur;
        Real cutoff_cur, noise_cur;

        sweep_table.file() >> maxm_cur >> minm_cur >> cutoff_cur >> niter_cur >> noise_cur;

        if (!sweep_table.file()) {
            // If we run out of values, use the last given values for the rest
            // of the sweeps.
            if (sw_phys < 1) {
                Error("Not enough sweeps provided");
            }

            for (; sw <= num_sweeps; sw++) {
                sweeps.setmaxm(sw, maxm_last);
                sweeps.setminm(sw, minm_last);
                sweeps.setcutoff(sw, cutoff_last);
                sweeps.setniter(sw, niter_last);
                sweeps.setnoise(sw, noise_last);
            }

            break;
        }

        maxm_last = maxm_cur;
        minm_last = minm_cur;
        cutoff_last = cutoff_cur;
        niter_last = niter_cur;
        noise_last = noise_cur;

        sw_phys++;

        if (sw_phys <= skip_sweeps) {
            // We haven't skipped enough yet, so ignore the last values.
            sw--;
        } else {
            sweeps.setmaxm(sw, maxm_last);
            sweeps.setminm(sw, minm_last);
            sweeps.setcutoff(sw, cutoff_last);
            sweeps.setniter(sw, niter_last);
            sweeps.setnoise(sw, noise_last);
        }
    }

    return sweeps;
}

void dmrg_sweep(IQMPS& psi, IQMPO const& H, InputGroup& sweep_table, int num_sweeps, int skip_sweeps, std::vector<IQMPS> ortho_wfs) {
    int N = psi.N();

    auto sweeps = make_sweeps(sweep_table, num_sweeps, skip_sweeps);
    println();
    println(sweeps);

    auto obs = LinRotObserver<IQTensor>(psi, H, num_sweeps);

    Real energy;
    if (ortho_wfs.empty()) {
        // Ground state.
        energy = dmrg(psi, H, sweeps, obs, "Quiet");
    } else {
        // Excited state.
        energy = dmrg(psi, H, ortho_wfs, sweeps, obs, {"Quiet", true, "Weight", 20.0});
    }

    println();
    printfln("E0 = %.15f", energy);
    printfln("dH2 = %.15f", obs.dH2());

    // Entanglement spectrum at middle bond.
    for (auto lambda : obs.middle_eigs()) {
        printfln("lambda = %.15f", lambda);
    }

    // Entanglement entropy.
    for (auto i : range1(N/2)) {
        printfln("SvN(%04d) = %.15f", i, obs.SvN(i));
    }
    for (auto i : range1(N/2)) {
        printfln("S2(%04d) = %.15f", i, obs.S2(i));
    }
    for (auto i : range1(N/2)) {
        printfln("Sinf(%04d) = %.15f", i, obs.Sinf(i));
    }

    for (auto i : range(ortho_wfs.size())) {
        auto x = overlapC(ortho_wfs[i], psi);
        printf("overlap(%02d) = %.15f", i, x.real());
        if (x.imag() > 1e-12) {
            printfln("+%.15fi", x.imag());
        } else if (x.imag() < -1e-12) {
            printfln("-%.15fi", -x.imag());
        } else {
            println();
        }
    }
}


void spatial_correlation(IQMPS& psi, std::string op1_name, std::string op2_name) {
    int N = psi.N();
    auto sites = psi.sites();

    Real *corr = new Real[N*N];

    for (auto i : range1(N)) {
        for (auto j : range1(i, N)) {
            corr[N*(i-1)+(j-1)] = 0.0;
        }
    }

    for (const auto& term_i : LinearRigidRotorSite::compound_op1(op1_name)) {
        for (const auto& term_j : LinearRigidRotorSite::compound_op1(op2_name)) {
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
                        printfln("WARNING: Complex correlation (%.15f)", zz.imag());
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
                        printfln("WARNING: Complex correlation (%.15f)", zz.imag());
                    }

                    corr[N*(i-1)+(j-1)] += zz.real();
                }
            }
        }
    }

    for (auto i : range1(N)) {
        for (auto j : range1(i, N)) {
            printfln("<%s%04d %s%04d> = %.15f", op1_name, i, op2_name, j, corr[N*(i-1)+(j-1)]);
        }
    }

    delete[] corr;
}

void spatial_correlation2(IQMPS& psi, std::string op1a_name, std::string op1b_name, std::string op2a_name, std::string op2b_name) {
    int N = psi.N();
    auto sites = psi.sites();

    Real *corr = new Real[N*N];

    for (auto i : range1(N)) {
        for (auto j : range1(i, N)) {
            corr[N*(i-1)+(j-1)] = 0.0;
        }
    }

    for (const auto& term_ia : LinearRigidRotorSite::compound_op1(op1a_name)) {
        for (const auto& term_ib : LinearRigidRotorSite::compound_op1(op1b_name)) {
            for (const auto& term_ja : LinearRigidRotorSite::compound_op1(op2a_name)) {
                for (const auto& term_jb : LinearRigidRotorSite::compound_op1(op2b_name)) {
                    for (auto i : range1(N)) {
                        psi.position(i);

                        auto op_ia = sites.op(term_ia.op, i);
                        auto op_ib = sites.op(term_ib.op, i);

                        // Diagonal.
                        {
                            auto op_i1a = op_ia;
                            auto op_i1b = op_ib;
                            auto op_i2a = sites.op(term_ja.op, i);
                            auto op_i2b = sites.op(term_jb.op, i);

                            auto C = psi.A(i) * op_i1a * prime(op_i1b) * prime(op_i2a, 2) * prime(op_i2b, 3) * dag(prime(psi.A(i), Site, 4));
                            auto zz = term_ia.k * term_ib.k * term_ja.k * term_jb.k * C.cplx();

                            if (std::abs(zz.imag()) > 1e-12) {
                                printfln("WARNING: Complex correlation (%.15f)", zz.imag());
                            }

                            corr[(N+1)*(i-1)] += zz.real();
                        }

                        // Off-diagonal.
                        auto idx = commonIndex(psi.A(i), psi.A(i+1), Link);
                        auto C = psi.A(i) * op_ia * prime(op_ib) * dag(prime(prime(psi.A(i), Site, 2), idx, 1));

                        for (auto j : range1(i+1, N)) {
                            auto op_ja = sites.op(term_ja.op, j);
                            auto op_jb = sites.op(term_jb.op, j);

                            if (j-1 > i) {
                                C *= psi.A(j-1);
                                C *= dag(prime(psi.A(j-1), Link));
                            }

                            auto D = C * psi.A(j);
                            D *= op_ja * prime(op_jb);
                            auto idx = commonIndex(psi.A(j), psi.A(j-1), Link);
                            D *= dag(prime(prime(psi.A(j), Site, 2), idx, 1));
                            auto zz = term_ia.k * term_ib.k * term_ja.k * term_jb.k * D.cplx();

                            if (std::abs(zz.imag()) > 1e-12) {
                                printfln("WARNING: Complex correlation (%.15f)", zz.imag());
                            }

                            corr[N*(i-1)+(j-1)] += zz.real();
                        }
                    }
                }
            }
        }
    }

    for (auto i : range1(N)) {
        for (auto j : range1(i, N)) {
            printfln("<%s%04d %s%04d %s%04d %s%04d> = %.15f", op1a_name, i, op1b_name, i, op2a_name, j, op2b_name, j, corr[N*(i-1)+(j-1)]);
        }
    }

    delete[] corr;
}

void run_analysis(IQMPS& psi) {
    int N = psi.N();
    auto sites = psi.sites();

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
    spatial_correlation(psi, "x", "x");
    spatial_correlation(psi, "z", "z");
    spatial_correlation2(psi, "z", "z", "z", "z");
}


void dump_coefficients(int l_max, IQMPS const& psi) {
    int N = psi.N();
    auto sites = psi.sites();

    // Start with {0, 1, 1, ..., 1, 1}.
    std::vector<int> conf(N);
    for (auto i : range1(2, N)) {
        conf[i-1]++;
    }

    for (;;) {
        bool done = false;
        for (auto i : range1(N)) {
            conf[i-1]++;
            if (conf[i-1] <= sites(i).m()) {
                break;
            } else {
                if (i == N) {
                    done = true;
                } else {
                    conf[i-1] = 1;
                }
            }
        }
        if (done) break;

        auto tensor = dag(setElt(sites(1)(conf[0]))) * psi.A(1);
        for (auto i : range1(2, N)) {
            tensor *= dag(setElt(sites(i)(conf[i-1]))) * psi.A(i);
        }
        auto coef = tensor.cplx();
        auto prob = std::norm(coef);

        // Only output configurations that have sufficient probability.
        if (prob < 1e-12) continue;

        print("coefficient: ");
        for (auto n : conf) {
            print(LinearRigidRotorSite::state_label(l_max, n), " ");
        }
        printf("%.15f", coef.real());
        if (coef.imag() > 1e-12) {
            printfln("+%.15fi", coef.imag());
        } else if (coef.imag() < -1e-12) {
            printfln("-%.15fi", -coef.imag());
        } else {
            println();
        }
    }
}


void run_sampling(int l_max, IQMPS& psi, int num_samples) {
    int N = psi.N();
    auto sites = psi.sites();

    std::random_device rd;
    std::mt19937 rng(rd());

    IQTensor cap;

    // Move the orthogonality center all the way to the first site so that the
    // other sites can be ignored when contracting.
    psi.position(1);

    for (auto sample_idx : range1(num_samples)) {
        std::vector<int> sample(N);

        for (auto i : range1(N)) {
            auto tensor = i == 1 ? psi.A(i) : cap * psi.A(i);
            auto rho = prime(dag(tensor), Site) * tensor;
            auto I = sites(i);
            auto Ip = prime(I);

            // Compute the weights of the outcomes. The first weight is a dummy
            // for indexing and should always remain zero.
            std::vector<Real> weights(I.m()+1);
            for (auto idx : range1(I.m())) {
                weights[idx] = rho.real(I(idx), Ip(idx));
            }

            // Sample.
            std::discrete_distribution<> d(weights.begin(), weights.end());
            sample[i-1] = d(rng);

            // Add site to cap.
            cap = dag(setElt(I(sample[i-1]))) * tensor;
        }

        print("sample ", sample_idx, ": ");
        for (auto n : sample) {
            print(LinearRigidRotorSite::state_label(l_max, n), " ");
        }
        println();
    }
}


IQIndex extract_link(IQTensor const& A, IQTensor const& B) {
    auto link = commonIndex(A, B, Link);
    if (link.dir() != dir(A, link)) {
        link.dag();
    }
    return link;
}

IQMPS embiggen(LinearRigidRotor const& sites1, IQMPS const& mps1, LinearRigidRotor const& sites2) {
    auto N = sites2.N();
    auto mps2 = IQMPS(sites2);

    {
        auto link_right = extract_link(mps1.A(1), mps1.A(2));
        auto A = IQTensor(sites2(1), link_right);
        for (auto idxR : range1(link_right.m())) {
            for (auto idxS1 : range1(sites1(1).m())) {
                auto val = mps1.A(1).cplx(sites1(1)(idxS1), link_right(idxR));
                if (val == Cplx(0, 0)) continue;
                auto label = LinearRigidRotorSite::state_label(sites1.l_max(), idxS1);
                A.set(sites2(1, label), link_right(idxR), val);
            }
        }
        mps2.setA(1, A);
    }

    for (auto i : range1(2, N-1)) {
        auto link_left = extract_link(mps1.A(i), mps1.A(i-1));
        auto link_right = extract_link(mps1.A(i), mps1.A(i+1));
        auto A = IQTensor(link_left, sites2(i), link_right);
        for (auto idxL : range1(link_left.m())) {
            for (auto idxR : range1(link_right.m())) {
                for (auto idxS1 : range1(sites1(i).m())) {
                    auto val = mps1.A(i).cplx(link_left(idxL), sites1(i)(idxS1), link_right(idxR));
                    if (val == Cplx(0, 0)) continue;
                    auto label = LinearRigidRotorSite::state_label(sites1.l_max(), idxS1);
                    A.set(link_left(idxL), sites2(i, label), link_right(idxR), val);
                }
            }
        }
        mps2.setA(i, A);
    }

    {
        auto link_left = extract_link(mps1.A(N), mps1.A(N-1));
        auto A = IQTensor(link_left, sites2(N));
        for (auto idxL : range1(link_left.m())) {
            for (auto idxS1 : range1(sites1(N).m())) {
                auto val = mps1.A(N).cplx(link_left(idxL), sites1(N)(idxS1));
                if (val == Cplx(0, 0)) continue;
                auto label = LinearRigidRotorSite::state_label(sites1.l_max(), idxS1);
                A.set(link_left(idxL), sites2(N, label), val);
            }
        }
        mps2.setA(N, A);
    }

    return mps2;
}
