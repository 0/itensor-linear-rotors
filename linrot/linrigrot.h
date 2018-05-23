#ifndef __LINROT_LINRIGROT_H
#define __LINROT_LINRIGROT_H

#include "itensor/mps/siteset.h"
#include "operators.h"


class LinearRigidRotorSite;

using LinearRigidRotor = itensor::BasicSiteSet<LinearRigidRotorSite>;


class LinearRigidRotorSite {
    int l_max;
    bool lp_sym;
    bool m_sym;
    itensor::IQIndex s;

protected:
    int m_offset(int l, int m) const {
        int m_offset = 0;
        for (int mp = -l_max; mp < m; mp++) {
            m_offset += (l_max-l%2)/2+1 - (std::abs(mp)+1-l%2)/2;
        }
        return m_offset;
    }

    int basis_idx(int l, int m) const {
        // The basis states are labelled by lp, m, and l quantum numbers, where
        // lp = l % 2 is the parity of l. The states are ordered so that l
        // changes most rapidly (in the sense of an innermost loop), m is next,
        // and lp changes most slowly (in the sense of an outer loop). Note
        // that while lp and m increment regularly, l decreases in steps of 2.
        // The order for lp and m is determined in the constructor when we
        // create the blocks for the IQIndex, and the order for l (within a
        // block) is determined here.
        //
        // For example, for l_max = 2, the states are ordered as follows:
        //     idx  lp   m  l
        //     1    0   -2  2
        //     2    0   -1  2
        //     3    0    0  2
        //     4    0    0  0
        //     5    0   +1  2
        //     6    0   +2  2
        //     7    1   -1  1
        //     8    1    0  1
        //     9    1   +1  1

        if (lp_sym) {
            if (m_sym) {
                return offset(s, itensor::findByQN(s, itensor::QN(l%2, m))) + (l_max-l)/2 + 1;
            } else {
                return offset(s, itensor::findByQN(s, itensor::QN(l%2))) + m_offset(l, m) + (l_max-l)/2 + 1;
            }
        } else {
            return (l%2)*(l_max/2+1)*(l_max-l_max%2+1) + m_offset(l, m) + (l_max-l)/2 + 1;
        }
    }

public:
    LinearRigidRotorSite(int n, itensor::Args const& args = itensor::Args::global()) {
        l_max = args.getInt("l_max");
        lp_sym = args.getBool("lp_sym", true);
        m_sym = args.getBool("m_sym", lp_sym);
        if (!lp_sym && m_sym) {
            itensor::Error("m_sym not allowed without lp_sym");
        }

        std::vector<itensor::IndexQN> iq;

        // The quantum numbers are total l parity (modulus 2) and total m
        // (modulus 1, meaning regular addition).
        if (lp_sym) {
            for (int lp = 0; lp <= 1; lp++) {
                if (m_sym) {
                    for (int m = -l_max; m <= l_max; m++) {
                        int size = (l_max - std::abs(m) + (std::abs(m) % 2 == lp) + (l_max % 2 == lp)) / 2;

                        if (size == 0) continue;

                        iq.push_back(itensor::IndexQN(itensor::Index(itensor::format("lp%dm%d:site%d", lp, m, n), size, itensor::Site),
                                     itensor::QN({lp, 2}, {m, 1})));
                    }
                } else {
                    int size = 2*((l_max+lp)/2)*((l_max+lp)/2+1) + (2*lp+1)*((l_max+lp)/2+1);

                    if (size == 0) continue;

                    iq.push_back(itensor::IndexQN(itensor::Index(itensor::format("lp%d:site%d", lp, n), size, itensor::Site),
                                 itensor::QN({lp, 2})));
                }
            }
        } else {
            int size = (l_max+1)*(l_max+1);

            iq.push_back(itensor::IndexQN(itensor::Index(itensor::format(":site%d", n), size, itensor::Site),
                         itensor::QN()));
        }

        s = itensor::IQIndex{itensor::nameint("rotor site=", n), std::move(iq)};
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

        if (opname == "B0") {
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
            Op.add_term(-2.0, "B0", "B0");
            Op.add_term(0.5, "B-", "B+");
            Op.add_term(0.5, "B+", "B-");
        } else if (opname == "D-D") {
            auto rx = args.getReal("rx");
            auto ry = args.getReal("ry");
            auto rz = args.getReal("rz");
            auto rz2 = rz*rz;
            auto rperp = itensor::Cplx(rx, ry);
            auto rperpc = std::conj(rperp);

            Op.add_term(1-3*rz2, "B0", "B0");
            Op.add_term(-0.25*(1-3*rz2), "B-", "B+");
            Op.add_term(-0.25*(1-3*rz2), "B+", "B-");
            Op.add_term(-0.75*rperp*rperp, "B-", "B-");
            Op.add_term(-0.75*rperpc*rperpc, "B+", "B+");
            Op.add_term(-1.5*rperp*rz, "B-", "B0");
            Op.add_term(-1.5*rperpc*rz, "B+", "B0");
            Op.add_term(-1.5*rperp*rz, "B0", "B-");
            Op.add_term(-1.5*rperpc*rz, "B0", "B+");
        } else {
            itensor::Error("Operator \"" + opname + "\" name not recognized");
        }

        return Op;
    }
};

#endif // __LINROT_LINRIGROT_H
