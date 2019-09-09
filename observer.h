#ifndef __LINROT_OBSERVER_H
#define __LINROT_OBSERVER_H

#include "itensor/mps/mps.h"
#include "itensor/mps/DMRGObserver.h"
#include "itensor/util/print.h"


template<class Tensor>
class LinRotObserver : public itensor::DMRGObserver<Tensor> {
    itensor::MPOt<Tensor> H;

    int nsweep;

    Tensor H2_contraction;
    itensor::Real dH2_;

    int max_num_eigs;

    std::vector<itensor::Real> middle_eigs_;
    std::vector<itensor::Real> SvNs;
    std::vector<itensor::Real> S2s;
    std::vector<itensor::Real> Sinfs;

public:
    LinRotObserver(itensor::MPSt<Tensor> const& psi,
                   itensor::MPOt<Tensor> const& H, int nsweep,
                   itensor::Args const& args = itensor::Args::global())
            : itensor::DMRGObserver<Tensor>(psi, args),
              H(H),
              nsweep(nsweep),
              dH2_(std::numeric_limits<itensor::Real>::infinity()),
              max_num_eigs(-1),
              middle_eigs_(16),
              SvNs(psi.N()-1),
              S2s(psi.N()-1),
              Sinfs(psi.N()-1)
            { }

    itensor::Real dH2() {
        return dH2_;
    }

    std::vector<itensor::Real> middle_eigs() {
        return middle_eigs_;
    }

    itensor::Real SvN(int i) {
        return SvNs[i-1];
    }

    itensor::Real S2(int i) {
        return S2s[i-1];
    }

    itensor::Real Sinf(int i) {
        return Sinfs[i-1];
    }

    void measure(itensor::Args const& args = itensor::Args::global()) {
        itensor::DMRGObserver<Tensor>::measure(args);

        auto psi = itensor::DMRGObserver<Tensor>::psi();
        auto N = psi.N();
        auto spectrum = itensor::DMRGObserver<Tensor>::spectrum();

        auto sw = args.getInt("Sweep");
        auto b = args.getInt("AtBond");
        auto ha = args.getInt("HalfSweep");
        auto energy = args.getReal("Energy");

        // First half-sweep.
        if (ha == 1) return;

        if (b == N/2) {
            // Middle of second half-sweep.
            std::vector<itensor::Real>::size_type idx = 0;
            for (auto p : spectrum.eigsKept()) {
                middle_eigs_[idx] = p;
                if (++idx >= middle_eigs_.size()) break;
            }
            for (; idx < middle_eigs_.size(); idx++) {
                middle_eigs_[idx] = 0.0;
            }
        }

        itensor::Real SvN_ = 0.0;
        itensor::Real S2_ = 0.0;
        itensor::Real max_eig = -1.0;
        for (auto p : spectrum.eigsKept()) {
            if (p < 1e-12) continue;

            SvN_ += -p*log(p);
            S2_ += p*p;

            if (p > max_eig) {
                max_eig = p;
            }
        }
        SvNs[b-1] = SvN_;
        S2s[b-1] = -log(S2_);
        Sinfs[b-1] = -log(max_eig);

        if (b+1 == N) {
            // Start of second half-sweep.
            max_num_eigs = spectrum.numEigsKept();

            H2_contraction = psi.A(b+1) * H.A(b+1) * prime(H.A(b+1)) * dag(prime(psi.A(b+1), 2));
        } else {
            // Rest of second half-sweep (including end).
            max_num_eigs = std::max(max_num_eigs, spectrum.numEigsKept());

            H2_contraction *= psi.A(b+1);
            H2_contraction *= H.A(b+1);
            H2_contraction *= prime(H.A(b+1));
            H2_contraction *= dag(prime(psi.A(b+1), 2));
        }

        if (b == 1) {
            // End of second half-sweep.
            H2_contraction *= psi.A(1);
            H2_contraction *= H.A(1);
            H2_contraction *= prime(H.A(1));
            H2_contraction *= dag(prime(psi.A(1), 2));

            if (std::abs(H2_contraction.cplx().imag()) > 1e-12) {
                itensor::println("WARNING: Complex H2!");
            } else if (H2_contraction.real() < 0.0) {
                itensor::println("WARNING: Negative H2!");
            } else if (H2_contraction.real() < energy*energy) {
                itensor::println("WARNING: H2 less than E^2!");
            }

            dH2_ = N*(H2_contraction.real()/(energy*energy) - 1.0);

            itensor::printfln("    dH2 after sweep %d is %.12f", sw, dH2_);
        }
    }

    bool checkDone(itensor::Args const& args = itensor::Args::global()) {
        auto sw = args.getInt("Sweep");
        auto maxm = args.getInt("Maxm");

        bool done = itensor::DMRGObserver<Tensor>::checkDone(args);
        done |= sw >= nsweep;

        if (done && max_num_eigs >= maxm) {
            itensor::println("WARNING: maxm reached on final sweep!");
        }

        return done;
    }
};

#endif // __LINROT_OBSERVER_H
