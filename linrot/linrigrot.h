#ifndef __LINROT_LINRIGROT_H
#define __LINROT_LINRIGROT_H

#include <map>

#include "itensor/mps/siteset.h"
#include "itensor/util/error.h"
#include "operators.h"


class LinearRigidRotorSite {
    static std::map<int,std::map<int,std::string> > state_map;

    std::map<std::tuple<int,int>,int> idx_map;

    int l_max;
    bool lp_sym;
    bool m_sym;
    itensor::IQIndex s;

protected:
    void set_args(itensor::Args const& args) {
        l_max = args.getInt("l_max");
        lp_sym = args.getBool("lp_sym", true);
        m_sym = args.getBool("m_sym", lp_sym);
        if (!lp_sym && m_sym) {
            itensor::Error("m_sym not allowed without lp_sym");
        }
    }

    int basis_idx(int l, int m) const {
        return idx_map.at({l, m});
    }

    std::map<std::tuple<int,int>,std::vector<std::tuple<int,int> > > blocks() const {
        // Mapping from quantum numbers identifying a symmetry block to the
        // states contained in that block.
        std::map<std::tuple<int,int>,std::vector<std::tuple<int,int> > > blocks;

        // The quantum numbers are l parity (modulus 2) and m (modulus 1,
        // meaning regular addition). Populate the blocks by enumerating all
        // possible one-rotor states.
        for (int l = 0; l <= l_max; l++) {
            int lp = l % 2;

            for (int m = -l; m <= l; m++) {
                if (lp_sym) {
                    if (m_sym) {
                        blocks[{lp, m}].push_back({l, m});
                    } else {
                        blocks[{lp, 0}].push_back({l, m});
                    }
                } else {
                    blocks[{0, 0}].push_back({l, m});
                }
            }
        }

        return blocks;
    }

    void populate_idx_map() {
        int state_idx = 1;

        for (auto block : blocks()) {
            for (auto state : block.second) {
                idx_map[state] = state_idx;
                state_idx++;
            }
        }
    }

    itensor::IQIndex construct_index(int n) const {
        std::vector<itensor::IndexQN> iq;

        for (auto block : blocks()) {
            int lp = std::get<0>(block.first);
            int m = std::get<1>(block.first);
            int size = block.second.size();

            if (lp_sym) {
                if (m_sym) {
                    auto index = itensor::Index(itensor::format("lp%dm%d:site%d", lp, m, n),
                                                size, itensor::Site);
                    auto qn = itensor::QN({lp, 2}, {m, 1});
                    iq.push_back(itensor::IndexQN(index, qn));
                } else {
                    auto index = itensor::Index(itensor::format("lp%d:site%d", lp, n),
                                                size, itensor::Site);
                    auto qn = itensor::QN({lp, 2});
                    iq.push_back(itensor::IndexQN(index, qn));
                }
            } else {
                auto index = itensor::Index(itensor::format(":site%d", n),
                                            size, itensor::Site);
                auto qn = itensor::QN();
                iq.push_back(itensor::IndexQN(index, qn));
            }
        }

        return itensor::IQIndex{itensor::nameint("rotor site=", n), std::move(iq)};
    }

    void populate_state_map() {
        if (state_map.count(l_max) == 0) {
            for (int l = 0; l <= l_max; l++) {
                for (int m = -l; m <= l; m++) {
                    std::ostringstream oss;
                    oss << 'l' << l << 'm' << m;
                    state_map[l_max][basis_idx(l, m)] = oss.str();
                }
            }
        }
    }

public:
    LinearRigidRotorSite(itensor::IQIndex s, itensor::Args const& args = itensor::Args::global()) : s(s) {
        set_args(args);
        populate_idx_map();
        populate_state_map();
    }

    LinearRigidRotorSite(int n, itensor::Args const& args = itensor::Args::global()) {
        set_args(args);
        s = construct_index(n);
        populate_idx_map();
        populate_state_map();
    }

    itensor::IQIndex index() const {
        return s;
    }

    itensor::IQIndexVal state(std::string const& state) {
        std::istringstream iss(state);
        char x1, x2;
        int l, m;

        iss >> x1 >> l >> x2 >> m;

        if (iss.fail() != 0 || x1 != 'l' || x2 != 'm' || l < 0 || l_max < l || m < -l || l < m) {
            itensor::Error("State " + state + " not recognized");

            return itensor::IQIndexVal{};
        }

        return s(basis_idx(l, m));
    }

    itensor::IQTensor op(std::string const& opname, itensor::Args const& args) const {
        auto sP = itensor::prime(s);
        auto Op = itensor::IQTensor(dag(s), sP);

        if (opname == "1") {
            for (int l = 0; l <= l_max; l++) {
                for (int m = -l; m <= l; m++) {
                    Op.set(s(basis_idx(l, m)),
                           sP(basis_idx(l, m)),
                           1);
                }
            }
        } else if (opname == "L2") {
            for (int l = 0; l <= l_max; l++) {
                for (int m = -l; m <= l; m++) {
                    Op.set(s(basis_idx(l, m)),
                           sP(basis_idx(l, m)),
                           l*(l+1));
                }
            }
        } else if (opname == "L+") {
            for (int l = 0; l <= l_max-1; l++) {
                for (int m = -l; m <= l; m++) {
                    Op.set(s(basis_idx(l, m)),
                           sP(basis_idx(l+1, m)),
                           std::sqrt(itensor::Real(l-m+1)*(l+m+1)/((2*l+1)*(2*l+3))));
                }
            }
        } else if (opname == "L-") {
            for (int l = 1; l <= l_max; l++) {
                for (int m = -l+1; m <= l-1; m++) {
                    Op.set(s(basis_idx(l, m)),
                           sP(basis_idx(l-1, m)),
                           std::sqrt(itensor::Real(l-m)*(l+m)/((2*l-1)*(2*l+1))));
                }
            }
        } else if (opname == "L+M+") {
            for (int l = 0; l <= l_max-1; l++) {
                for (int m = -l; m <= l; m++) {
                    Op.set(s(basis_idx(l, m)),
                           sP(basis_idx(l+1, m+1)),
                           -std::sqrt(itensor::Real(l+m+1)*(l+m+2)/((2*l+1)*(2*l+3))));
                }
            }
        } else if (opname == "L+M-") {
            for (int l = 0; l <= l_max-1; l++) {
                for (int m = -l; m <= l; m++) {
                    Op.set(s(basis_idx(l, m)),
                           sP(basis_idx(l+1, m-1)),
                           -std::sqrt(itensor::Real(l-m+1)*(l-m+2)/((2*l+1)*(2*l+3))));
                }
            }
        } else if (opname == "L-M+") {
            for (int l = 1; l <= l_max; l++) {
                for (int m = -l; m <= l-2; m++) {
                    Op.set(s(basis_idx(l, m)),
                           sP(basis_idx(l-1, m+1)),
                           std::sqrt(itensor::Real(l-m-1)*(l-m)/((2*l-1)*(2*l+1))));
                }
            }
        } else if (opname == "L-M-") {
            for (int l = 1; l <= l_max; l++) {
                for (int m = -l+2; m <= l; m++) {
                    Op.set(s(basis_idx(l, m)),
                           sP(basis_idx(l-1, m-1)),
                           std::sqrt(itensor::Real(l+m-1)*(l+m)/((2*l-1)*(2*l+1))));
                }
            }
        } else {
            itensor::Error("Operator \"" + opname + "\" name not recognized");
        }

        return Op;
    }

    static OneSiteOperator compound_op1(std::string const& opname) {
        OneSiteOperator Op;

        if (opname == "B0" || opname == "z") {
            Op.add_term(1.0, "L+");
            Op.add_term(1.0, "L-");
        } else if (opname == "B+") {
            Op.add_term(1.0, "L+M+");
            Op.add_term(1.0, "L-M+");
        } else if (opname == "B-") {
            Op.add_term(-1.0, "L+M-");
            Op.add_term(-1.0, "L-M-");
        } else if (opname == "x") {
            Op.add_term(0.5, "L+M+");
            Op.add_term(0.5, "L-M+");
            Op.add_term(-0.5, "L+M-");
            Op.add_term(-0.5, "L-M-");
        } else {
            // Pass unrecognized operators through to the real op method.
            Op.add_term(1.0, opname);
        }

        return Op;
    }

    static TwoSiteOperator compound_op2(std::string const& opname, itensor::Args const& args = itensor::Args::global()) {
        TwoSiteOperator Op(&compound_op1);

        if (opname == "dot product") {
            Op.add_term(1.0, "B0", "B0");
            Op.add_term(0.5, "B-", "B+");
            Op.add_term(0.5, "B+", "B-");
        } else if (opname == "D-D lin") {
            auto a = args.getReal("a", -3.0);
            Op.add_term(1.0+a, "B0", "B0");
            Op.add_term(0.5, "B-", "B+");
            Op.add_term(0.5, "B+", "B-");
        } else if (opname == "D-D") {
            auto a = args.getReal("a", -3.0);
            auto rx = args.getReal("rx");
            auto ry = args.getReal("ry");
            auto rz = args.getReal("rz");
            auto rz2 = rz*rz;
            auto rperp = itensor::Cplx(rx, ry);
            auto rperpc = std::conj(rperp);
            Op.add_term(1.0+a*rz2, "B0", "B0");
            Op.add_term(0.25*(2.0+a-a*rz2), "B-", "B+");
            Op.add_term(0.25*(2.0+a-a*rz2), "B+", "B-");
            Op.add_term(0.25*a*rperp*rperp, "B-", "B-");
            Op.add_term(0.25*a*rperpc*rperpc, "B+", "B+");
            Op.add_term(0.5*a*rperp*rz, "B-", "B0");
            Op.add_term(0.5*a*rperpc*rz, "B+", "B0");
            Op.add_term(0.5*a*rperp*rz, "B0", "B-");
            Op.add_term(0.5*a*rperpc*rz, "B0", "B+");
        } else {
            itensor::Error("Operator \"" + opname + "\" name not recognized");
        }

        return Op;
    }

    static std::string state_label(int l_max, int idx) {
        return state_map.at(l_max).at(idx);
    }
};


class LinearRigidRotor : public itensor::SiteSet {
    int l_max_;
    bool lp_sym_;
    bool m_sym_;

public:
    LinearRigidRotor() { }

    LinearRigidRotor(int N, itensor::Args const& args = itensor::Args::global()) {
        l_max_ = args.getInt("l_max");
        lp_sym_ = args.getBool("lp_sym", true);
        m_sym_ = args.getBool("m_sym", lp_sym_);

        auto sites = itensor::SiteStore(N);
        for (auto i : itensor::range1(N)) {
            sites.set(i, LinearRigidRotorSite(i, args));
        }
        SiteSet::init(std::move(sites));
    }

    int l_max() const {
        return l_max_;
    }

    bool lp_sym() const {
        return lp_sym_;
    }

    bool m_sym() const {
        return m_sym_;
    }

    void write(std::ostream& s) const {
        if (N() == 0) {
            itensor::Error("Refusing to write empty SiteSet");
        }

        itensor::write(s, l_max_);
        itensor::write(s, lp_sym_);
        itensor::write(s, m_sym_);
        itensor::write(s, N());

        for (auto i : itensor::range1(N())) {
            this->operator()(i).write(s);
        }
    }

    void read(std::istream& s) {
        l_max_ = itensor::read<int>(s);
        lp_sym_ = itensor::read<bool>(s);
        m_sym_ = itensor::read<bool>(s);
        int N = itensor::read<int>(s);

        if (N == 0) {
            itensor::Error("Refusing to read empty SiteSet");
        }

        auto args = itensor::Args{"l_max", l_max_, "lp_sym", lp_sym_, "m_sym", m_sym_};
        auto sites = itensor::SiteStore(N);
        for (auto i : itensor::range1(N)) {
            auto I = itensor::IQIndex{};
            I.read(s);
            sites.set(i, LinearRigidRotorSite(I, args));
        }
        SiteSet::init(std::move(sites));
    }
};

#endif // __LINROT_LINRIGROT_H
