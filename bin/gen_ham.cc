#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "linrot/operators.h"
#include "argparse.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s [--pbc] [--geom <G>]"
                          " [--field <F> [--field-linear]] [--anisotropy <A>]"
                          " {-R <R> | -g <G>} [--sociability <S>]"
                          " --mpo-cutoff <C> --sites <S>"
                          " --ham-out <H>", argv[0]);
        return 1;
    }

    ArgumentParser parser;
    parser.add("--pbc", ArgType::Flag, {"required", false});
    parser.add("--geom", ArgType::String, {"required", false});
    parser.add("--field", ArgType::Real, {"required", false});
    parser.add("--field-linear", ArgType::Flag, {"required", false});
    parser.add("--anisotropy", ArgType::Real, {"required", false});
    parser.add("-R", ArgType::Real, {"required", false});
    parser.add("-g", ArgType::Real, {"required", false});
    parser.add("--sociability", ArgType::Int, {"required", false});
    parser.add("--mpo-cutoff", ArgType::Real);
    parser.add("--sites", ArgType::String);
    parser.add("--ham-out", ArgType::String);
    auto args = parser.parse(argc, argv);

    Real inter;
    if (!args.defined("R") && !args.defined("g")) {
        println("One of -R or -g must be provided");
        return 1;
    } else if (args.defined("R") && args.defined("g")) {
        println("Only one of -R or -g may be provided");
        return 1;
    } else if (args.defined("R")) {
        Real R = args.getReal("R");
        inter = 1.0/(R*R*R);
    } else {
        Real g = args.getReal("g");
        inter = g;
    }

    bool pbc = args.getBool("pbc", false);
    auto geom_in_path = args.getString("geom", "");
    Real field = args.getReal("field", 0.0);
    bool field_linear = args.getBool("field-linear", false);
    Real anisotropy = args.getReal("anisotropy", -3.0);
    int sociability = args.getInt("sociability", -1);
    Real mpo_cutoff = args.getReal("mpo-cutoff");
    auto sites_in_path = args.getString("sites");
    auto H_out_path = args.getString("ham-out");

    auto sites = readFromFile<LinearRigidRotor>(sites_in_path);
    auto N = length(sites);

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
        if (sociability >= 0) {
            println("No support for --sociability with --geom-in-path");
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

    cpu_time time;
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
                Real k = inter/(r*r*r);
                add_operator(ampo, LinearRigidRotorSite::compound_op2("D-D", {"a", anisotropy, "rx", rx/r, "ry", ry/r, "rz", rz/r}), i, j, k);
            } else {
                int d = j-i;
                if (pbc && j-i > N/2) {
                    // Minimum image convention for periodic boundary conditions.
                    d = N - d;
                }
                if (sociability >= 0 && d > sociability) continue;
                Real k = inter/(d*d*d);
                add_operator(ampo, LinearRigidRotorSite::compound_op2("D-D lin", {"a", anisotropy}), i, j, k);
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
    auto H = toMPO(ampo, {"Cutoff", mpo_cutoff});
    printfln("cputime = %.15e", time.sincemark().time);
    printfln("walltime = %.15e", time.sincemark().wall);
    printfln("mean MPO bond dimension: %.15e", averageLinkDim(H));

    if (geom_nonlinear) {
        delete[] geom;
    }

    writeToFile(H_out_path, H);

    return 0;
}
