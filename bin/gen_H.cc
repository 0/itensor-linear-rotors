#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "linrot/operators.h"
#include "argparse.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s [--pbc] [--geom-in-path <G>]"
                          " [--field <F> [--field-linear]] -R <R>"
                          " --mpo-cutoff <C> --sites-in-path <S>"
                          " --H-out-path <H>", argv[0]);
        return 1;
    }

    ArgumentParser parser;
    parser.add("--pbc", ArgType::Flag, {"required", false});
    parser.add("--geom-in-path", ArgType::String, {"required", false});
    parser.add("--field", ArgType::Real, {"required", false});
    parser.add("--field-linear", ArgType::Flag, {"required", false});
    parser.add("-R", ArgType::Real);
    parser.add("--mpo-cutoff", ArgType::Real);
    parser.add("--sites-in-path", ArgType::String);
    parser.add("--H-out-path", ArgType::String);
    auto args = parser.parse(argc, argv);

    bool pbc = args.getBool("pbc", false);
    auto geom_in_path = args.getString("geom-in-path", "");
    Real field = args.getReal("field", 0.0);
    bool field_linear = args.getBool("field-linear", false);
    Real R = args.getReal("R");
    Real mpo_cutoff = args.getReal("mpo-cutoff");
    auto sites_in_path = args.getString("sites-in-path");
    auto H_out_path = args.getString("H-out-path");

    auto sites = readFromFile<LinearRigidRotor>(sites_in_path);
    auto N = sites.N();

    bool geom_nonlinear = false;
    double *geom;
    if (!geom_in_path.empty()) {
        geom_nonlinear = true;
        geom = new double[3*N];

        std::ifstream geom_ifs(geom_in_path);
        for (auto i : range(N)) {
            geom_ifs >> geom[3*i+0] >> geom[3*i+1] >> geom[3*i+2];
        }
    }

    if (geom_nonlinear) {
        if (pbc) {
            println("No support for --pbc with --geom-in-path");
            return 1;
        }
        if (field_linear) {
            println("No support for --field-linear with --geom-in-path");
            return 1;
        }
        if (sites.m_sym()) {
            println("Sites must not have m symmetry with --geom-in-path");
            return 1;
        }
    }
    if (field != 0.0) {
        if (sites.lp_sym()) {
            println("Sites must not have lp symmetry with --field");
            return 1;
        }
    }

    auto ampo = AutoMPO(sites);
    // Rotational energy.
    for (auto i : range1(N)) {
        ampo += "L2", i;
    }
    // Potential energy.
    for (auto i : range1(N)) {
        for (auto j : range1(i+1, N)) {
            if (geom_nonlinear) {
                Real rx = geom[3*(j-1)+0] - geom[3*(i-1)+0];
                Real ry = geom[3*(j-1)+1] - geom[3*(i-1)+1];
                Real rz = geom[3*(j-1)+2] - geom[3*(i-1)+2];
                Real r = std::sqrt(rx*rx + ry*ry + rz*rz);
                Real d = R*r;
                Real k = 1.0/(d*d*d);
                add_operator(ampo, LinearRigidRotorSite::compound_op2("D-D", {"rx", rx/r, "ry", ry/r, "rz", rz/r}), i, j, k);
            } else {
                Real d = R*(j-i);
                if (pbc && j-i > N/2) {
                    // Minimum image convention for periodic boundary conditions.
                    d = R*N - d;
                }
                Real k = 1.0/(d*d*d);
                add_operator(ampo, LinearRigidRotorSite::compound_op2("D-D lin"), i, j, k);
            }
        }
    }
    // Transverse field.
    if (field != 0.0) {
        for (auto i : range1(N)) {
            Real k = field;
            if (field_linear) {
                k *= (double)(N-i)/(N-1);
            }
            add_operator(ampo, LinearRigidRotorSite::compound_op1("x"), i, k);
        }
    }
    auto H = toMPO<IQTensor>(ampo, {"Cutoff=", mpo_cutoff});
    printfln("mean MPO bond dimension: %f", averageM(H));

    if (geom_nonlinear) {
        delete[] geom;
    }

    writeToFile(H_out_path, H);

    return 0;
}
