#include "observer.h"
#include "dmrg.h"


std::string cplx2str(Cplx x) {
    std::ostringstream oss;
    printf(oss, "%.15e", x.real());
    if (x.imag() > 1e-12) {
        printf(oss, "+%.15ei", x.imag());
    } else if (x.imag() < -1e-12) {
        printf(oss, "-%.15ei", -x.imag());
    }
    return oss.str();
}


Sweeps make_sweeps(InputGroup& sweep_table, int num_sweeps, int skip_sweeps, Real mps_cutoff) {
    if (!sweep_table.GotoGroup()) {
        Error("Table \"" + sweep_table.name() + "\" not found");
    }

    // Skip the header.
    sweep_table.SkipLine();

    auto sweeps = Sweeps(num_sweeps);

    int sw_phys = 0;
    int maxdim_last, mindim_last, niter_last;
    Real cutoff_last, noise_last;

    for (int sw = 1; sw <= num_sweeps; sw++) {
        int maxdim_cur, mindim_cur, niter_cur;
        Real cutoff_cur, noise_cur;

        sweep_table.file() >> maxdim_cur >> mindim_cur >> cutoff_cur >> niter_cur >> noise_cur;

        if (!sweep_table.file()) {
            // If we run out of values, use the last given values for the rest
            // of the sweeps.
            if (sw_phys < 1) {
                Error("Not enough sweeps provided");
            }

            for (; sw <= num_sweeps; sw++) {
                sweeps.setmaxdim(sw, maxdim_last);
                sweeps.setmindim(sw, mindim_last);
                sweeps.setcutoff(sw, cutoff_last);
                sweeps.setniter(sw, niter_last);
                sweeps.setnoise(sw, noise_last);
            }

            break;
        }

        maxdim_last = maxdim_cur;
        mindim_last = mindim_cur;
        cutoff_last = cutoff_cur;
        niter_last = niter_cur;
        noise_last = noise_cur;

        sw_phys++;

        if (sw_phys <= skip_sweeps) {
            // We haven't skipped enough yet, so ignore the last values.
            sw--;
        } else {
            sweeps.setmaxdim(sw, maxdim_last);
            sweeps.setmindim(sw, mindim_last);
            sweeps.setcutoff(sw, cutoff_last);
            sweeps.setniter(sw, niter_last);
            sweeps.setnoise(sw, noise_last);
        }
    }

    if (mps_cutoff >= 0.0) {
        sweeps.cutoff() = mps_cutoff;
    }

    return sweeps;
}

void dmrg_sweep(MPS& psi, MPO const& H, InputGroup& sweep_table, int num_sweeps, int skip_sweeps, std::vector<MPS> ortho_wfs, Real mps_cutoff) {
    auto psi0 = psi;
    int N = length(psi);

    auto sweeps = make_sweeps(sweep_table, num_sweeps, skip_sweeps, mps_cutoff);
    println();
    println(sweeps);

    auto obs = LinRotObserver(psi, H);

    Real energy;
    if (ortho_wfs.empty()) {
        // Ground state.
        energy = dmrg(psi, H, sweeps, obs, "Quiet");
    } else {
        // Excited state.
        energy = dmrg(psi, H, ortho_wfs, sweeps, obs, {"Quiet", true, "Weight", 20.0});
    }

    println();
    printfln("E0 = %.15e", energy);
    printfln("dH2 = %.15e", obs.dH2());

    // Entanglement spectrum at middle bond.
    for (auto lambda : obs.middle_eigs()) {
        printfln("lambda = %.15e", lambda);
    }

    // Entanglement entropy.
    for (auto i : range1(N/2)) {
        printfln("SvN(%04d) = %.15e", i, obs.SvN(i));
    }
    for (auto i : range1(N/2)) {
        printfln("S2(%04d) = %.15e", i, obs.S2(i));
    }
    for (auto i : range1(N/2)) {
        printfln("Sinf(%04d) = %.15e", i, obs.Sinf(i));
    }

    // Overlap with previous state.
    printfln("overlap0 = %s", cplx2str(innerC(psi0, psi)));

    // Overlaps with orthogonal states.
    for (auto i : range(ortho_wfs.size())) {
        auto x = innerC(ortho_wfs[i], psi);
        printfln("overlap(%02d) = %s", i, cplx2str(x));
    }
}


void spatial_correlation(LinearRigidRotor const& sites, MPS& psi, std::string op1_name, std::string op2_name) {
    int N = length(psi);

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

                auto op_i = op(sites, term_i.op, i);

                // Diagonal.
                {
                    auto op_i1 = op_i;
                    auto op_i2 = op(sites, term_j.op, i);

                    auto C = psi(i) * op_i1 * prime(op_i2) * dag(prime(psi(i), 2, "Site"));
                    auto zz = term_i.k * term_j.k * C.cplx();

                    if (std::abs(zz.imag()) > 1e-12) {
                        printfln("WARNING: Complex correlation (%.15e)", zz.imag());
                    }

                    corr[(N+1)*(i-1)] += zz.real();
                }

                // Off-diagonal.
                auto idx = commonIndex(psi(i), psi(i+1), "Link");
                auto C = psi(i) * op_i * dag(prime(prime(psi(i), "Site"), idx));

                for (auto j : range1(i+1, N)) {
                    auto op_j = op(sites, term_j.op, j);

                    if (j-1 > i) {
                        C *= psi(j-1);
                        C *= dag(prime(psi(j-1), "Link"));
                    }

                    auto D = C * psi(j);
                    D *= op_j;
                    auto idx = commonIndex(psi(j), psi(j-1), "Link");
                    D *= dag(prime(prime(psi(j), "Site"), idx));
                    auto zz = term_i.k * term_j.k * D.cplx();

                    if (std::abs(zz.imag()) > 1e-12) {
                        printfln("WARNING: Complex correlation (%.15e)", zz.imag());
                    }

                    corr[N*(i-1)+(j-1)] += zz.real();
                }
            }
        }
    }

    for (auto i : range1(N)) {
        for (auto j : range1(i, N)) {
            printfln("<%s%04d %s%04d> = %.15e", op1_name, i, op2_name, j, corr[N*(i-1)+(j-1)]);
        }
    }

    delete[] corr;
}

void spatial_correlation2(LinearRigidRotor const& sites, MPS& psi, std::string op1a_name, std::string op1b_name, std::string op2a_name, std::string op2b_name) {
    int N = length(psi);

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

                        auto op_ia = op(sites, term_ia.op, i);
                        auto op_ib = op(sites, term_ib.op, i);

                        // Diagonal.
                        {
                            auto op_i1a = op_ia;
                            auto op_i1b = op_ib;
                            auto op_i2a = op(sites, term_ja.op, i);
                            auto op_i2b = op(sites, term_jb.op, i);

                            auto C = psi(i) * op_i1a * prime(op_i1b) * prime(op_i2a, 2) * prime(op_i2b, 3) * dag(prime(psi(i), 4, "Site"));
                            auto zz = term_ia.k * term_ib.k * term_ja.k * term_jb.k * C.cplx();

                            if (std::abs(zz.imag()) > 1e-12) {
                                printfln("WARNING: Complex correlation (%.15e)", zz.imag());
                            }

                            corr[(N+1)*(i-1)] += zz.real();
                        }

                        // Off-diagonal.
                        auto idx = commonIndex(psi(i), psi(i+1), "Link");
                        auto C = psi(i) * op_ia * prime(op_ib) * dag(prime(prime(psi(i), 2, "Site"), idx));

                        for (auto j : range1(i+1, N)) {
                            auto op_ja = op(sites, term_ja.op, j);
                            auto op_jb = op(sites, term_jb.op, j);

                            if (j-1 > i) {
                                C *= psi(j-1);
                                C *= dag(prime(psi(j-1), "Link"));
                            }

                            auto D = C * psi(j);
                            D *= op_ja * prime(op_jb);
                            auto idx = commonIndex(psi(j), psi(j-1), "Link");
                            D *= dag(prime(prime(psi(j), 2, "Site"), idx));
                            auto zz = term_ia.k * term_ib.k * term_ja.k * term_jb.k * D.cplx();

                            if (std::abs(zz.imag()) > 1e-12) {
                                printfln("WARNING: Complex correlation (%.15e)", zz.imag());
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
            printfln("<%s%04d %s%04d %s%04d %s%04d> = %.15e", op1a_name, i, op1b_name, i, op2a_name, j, op2b_name, j, corr[N*(i-1)+(j-1)]);
        }
    }

    delete[] corr;
}

void run_analysis(LinearRigidRotor const& sites, MPS& psi) {
    int N = length(psi);

    // Orientational correlation.
    auto OC_ampo = AutoMPO(sites);
    for (auto i : range1(N)) {
        for (auto j : range1(i+1, N)) {
            add_operator(OC_ampo, LinearRigidRotorSite::compound_op2("dot product"), i, j, 1.0);
        }
    }
    auto OC_op = toMPO(OC_ampo, {"Exact", true});
    auto OC = innerC(psi, OC_op, psi)*2.0/(N*(N-1));
    println("orientational correlation = ", cplx2str(OC));

    // Spatial correlations.
    spatial_correlation(sites, psi, "x", "x");
    spatial_correlation(sites, psi, "z", "z");
    spatial_correlation2(sites, psi, "z", "z", "z", "z");
}


void dump_coefficients(LinearRigidRotor const& sites, MPS const& psi) {
    int l_max = sites.l_max();
    int N = length(psi);

    // Start with {0, 1, 1, ..., 1, 1}.
    std::vector<int> conf(N);
    for (auto i : range1(2, N)) {
        conf[i-1]++;
    }

    for (;;) {
        bool done = false;
        for (auto i : range1(N)) {
            conf[i-1]++;
            if (conf[i-1] <= dim(sites(i))) {
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

        auto tensor = dag(setElt(sites(1)(conf[0]))) * psi(1);
        for (auto i : range1(2, N)) {
            tensor *= dag(setElt(sites(i)(conf[i-1]))) * psi(i);
        }
        auto coef = tensor.cplx();
        auto prob = std::norm(coef);

        // Only output configurations that have sufficient probability.
        if (prob < 1e-12) continue;

        print("coefficient: ");
        for (auto n : conf) {
            print(LinearRigidRotorSite::state_label(l_max, n), " ");
        }
        println(cplx2str(coef));
    }
}


void run_sampling(LinearRigidRotor const& sites, MPS& psi, int num_samples) {
    int l_max = sites.l_max();
    int N = length(psi);

    std::random_device rd;
    std::mt19937 rng(rd());

    ITensor cap;

    // Move the orthogonality center all the way to the first site so that the
    // other sites can be ignored when contracting.
    psi.position(1);

    for (auto sample_idx : range1(num_samples)) {
        std::vector<int> sample(N);

        for (auto i : range1(N)) {
            auto tensor = i == 1 ? psi(i) : cap * psi(i);
            auto rho = prime(dag(tensor), "Site") * tensor;
            auto I = sites(i);
            auto Ip = prime(I);

            // Compute the weights of the outcomes. The first weight is a dummy
            // for indexing and should always remain zero.
            std::vector<Real> weights(dim(I)+1);
            for (auto idx : range1(dim(I))) {
                weights[idx] = rho.cplx(I(idx), Ip(idx)).real();
            }

            // Sample.
            std::discrete_distribution<> d(weights.begin(), weights.end());
            sample[i-1] = d(rng);

            // Add site to cap.
            cap = dag(setElt(I(sample[i-1]))) * tensor;
        }

        print("sample ", sample_idx, ":");
        for (auto n : sample) {
            print(" ", LinearRigidRotorSite::state_label(l_max, n));
        }
        println();
    }
}


MPS embiggen(LinearRigidRotor const& sites1, MPS const& mps1, LinearRigidRotor const& sites2) {
    auto N = length(sites2);
    auto mps2 = MPS(sites2);

    {
        auto link_right = rightLinkIndex(mps1, 1);
        auto A = ITensor(sites2(1), link_right);
        for (auto idxR : range1(dim(link_right))) {
            for (auto idxS1 : range1(dim(sites1(1)))) {
                auto val = mps1(1).cplx(sites1(1)(idxS1), link_right(idxR));
                if (val == Cplx(0, 0)) continue;
                auto label = LinearRigidRotorSite::state_label(sites1.l_max(), idxS1);
                A.set(sites2(1, label), link_right(idxR), val);
            }
        }
        mps2.set(1, A);
    }

    for (auto i : range1(2, N-1)) {
        auto link_left = leftLinkIndex(mps1, i);
        auto link_right = rightLinkIndex(mps1, i);
        auto A = ITensor(link_left, sites2(i), link_right);
        for (auto idxL : range1(dim(link_left))) {
            for (auto idxR : range1(dim(link_right))) {
                for (auto idxS1 : range1(dim(sites1(i)))) {
                    auto val = mps1(i).cplx(link_left(idxL), sites1(i)(idxS1), link_right(idxR));
                    if (val == Cplx(0, 0)) continue;
                    auto label = LinearRigidRotorSite::state_label(sites1.l_max(), idxS1);
                    A.set(link_left(idxL), sites2(i, label), link_right(idxR), val);
                }
            }
        }
        mps2.set(i, A);
    }

    {
        auto link_left = leftLinkIndex(mps1, N);
        auto A = ITensor(link_left, sites2(N));
        for (auto idxL : range1(dim(link_left))) {
            for (auto idxS1 : range1(dim(sites1(N)))) {
                auto val = mps1(N).cplx(link_left(idxL), sites1(N)(idxS1));
                if (val == Cplx(0, 0)) continue;
                auto label = LinearRigidRotorSite::state_label(sites1.l_max(), idxS1);
                A.set(link_left(idxL), sites2(N, label), val);
            }
        }
        mps2.set(N, A);
    }

    return mps2;
}
