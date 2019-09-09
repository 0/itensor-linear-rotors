#ifndef __LINROT_OBSERVER_H
#define __LINROT_OBSERVER_H

#include "itensor/mps/DMRGObserver.h"
#include "itensor/mps/mps.h"
#include "itensor/mps/mpo.h"
#include "itensor/util/print.h"


class LinRotObserver : public itensor::DMRGObserver {
    itensor::MPO H;

    itensor::ITensor H2_contraction;
    itensor::Real dH2_;

    int max_num_eigs;

    std::vector<itensor::Real> middle_eigs_;
    std::vector<itensor::Real> SvNs;
    std::vector<itensor::Real> S2s;
    std::vector<itensor::Real> Sinfs;

public:
    LinRotObserver(itensor::MPS const& psi, itensor::MPO const& H,
                   itensor::Args const& args = itensor::Args::global())
            : itensor::DMRGObserver(psi, args),
              H(H),
              dH2_(std::numeric_limits<itensor::Real>::infinity()),
              max_num_eigs(-1),
              middle_eigs_(16),
              SvNs(length(psi)-1),
              S2s(length(psi)-1),
              Sinfs(length(psi)-1)
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
        itensor::DMRGObserver::measure(args);

        auto psi = itensor::DMRGObserver::psi();
        auto N = length(psi);
        auto spectrum = itensor::DMRGObserver::spectrum();

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

            H2_contraction = psi(b+1) * H(b+1) * prime(H(b+1)) * dag(prime(psi(b+1), 2));
        } else {
            // Rest of second half-sweep (including end).
            max_num_eigs = std::max(max_num_eigs, spectrum.numEigsKept());

            H2_contraction *= psi(b+1);
            H2_contraction *= H(b+1);
            H2_contraction *= prime(H(b+1));
            H2_contraction *= dag(prime(psi(b+1), 2));
        }

        if (b == 1) {
            // End of second half-sweep.
            H2_contraction *= psi(1);
            H2_contraction *= H(1);
            H2_contraction *= prime(H(1));
            H2_contraction *= dag(prime(psi(1), 2));

            if (std::abs(H2_contraction.cplx().imag()) > 1e-12) {
                itensor::printfln("WARNING: Complex H2 (%.15e)", H2_contraction.cplx().imag());
            } else if (H2_contraction.cplx().real() < 0.0) {
                itensor::printfln("WARNING: Negative H2 (%.15e)", H2_contraction.cplx().real());
            } else if (H2_contraction.cplx().real() < energy*energy) {
                itensor::printfln("WARNING: H2 less than E^2 (%.15e < %.15e)", H2_contraction.cplx().real(), energy*energy);
            }

            dH2_ = N*(H2_contraction.cplx().real()/(energy*energy) - 1.0);

            itensor::printfln("    dH2 after sweep %d is %.15e", sw, dH2_);
        }
    }

    bool checkDone(itensor::Args const& args = itensor::Args::global()) {
        auto sw = args.getInt("Sweep");
        auto nsweep = args.getInt("NSweep");
        auto maxdim = args.getInt("MaxDim");

        bool done = itensor::DMRGObserver::checkDone(args);
        done |= sw >= nsweep;

        if (done && max_num_eigs >= maxdim) {
            itensor::println("WARNING: maxdim reached on final sweep");
        }

        return done;
    }
};

#endif // __LINROT_OBSERVER_H
