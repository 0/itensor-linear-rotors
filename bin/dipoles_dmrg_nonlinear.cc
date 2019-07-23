#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "linrot/operators.h"
#include "argparse.h"

#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s --geom <G> -R <R> -N <N> --l-max <L>"
                          " --sweep-table <T> --dH2-goal <G>"
                          " --sweeps-min <S> --sweeps-max <S>", argv[0]);
        return 0;
    }

    ArgumentParser parser;
    parser.add("--geom", ArgType::String);
    parser.add("-R", ArgType::Real);
    parser.add("-N", ArgType::Int);
    parser.add("--l-max", ArgType::Int);
    parser.add("--sweep-table", ArgType::String);
    parser.add("--dH2-goal", ArgType::Real);
    parser.add("--sweeps-min", ArgType::Int);
    parser.add("--sweeps-max", ArgType::Int);
    auto args = parser.parse(argc, argv);

    Real R = args.getReal("R");
    auto geom_path = args.getString("geom");
    int N = args.getInt("N");
    int l_max = args.getInt("l-max");
    auto sweep_table = InputGroup(args.getString("sweep-table"), "sweeps");
    Real dH2_goal = args.getReal("dH2-goal");
    int sweeps_min = args.getInt("sweeps-min");
    int sweeps_max = args.getInt("sweeps-max");

    double *geom = new double[3*N];

    {
        std::ifstream geom_ifs(geom_path);
        for (auto i : range(N)) {
            geom_ifs >> geom[3*i+0] >> geom[3*i+1] >> geom[3*i+2];
        }
    }

    auto sites = LinearRigidRotor(N, {"l_max", l_max, "m_sym", false});

    auto ampo = AutoMPO(sites);
    // Rotational energy.
    for (auto i : range1(N)) {
        ampo += "L2", i;
    }
    // Potential energy.
    for (auto i : range1(N)) {
        for (auto j : range1(i+1, N)) {
            Real rx = geom[3*(j-1)+0] - geom[3*(i-1)+0];
            Real ry = geom[3*(j-1)+1] - geom[3*(i-1)+1];
            Real rz = geom[3*(j-1)+2] - geom[3*(i-1)+2];
            Real r = std::sqrt(rx*rx + ry*ry + rz*rz);
            Real d = R*r;
            Real k = 1.0/(d*d*d);

            add_operator(ampo, LinearRigidRotorSite::compound_op2("D-D", {"rx", rx/r, "ry", ry/r, "rz", rz/r}), i, j, k);
        }
    }
    auto H = toMPO<IQTensor>(ampo, {"Cutoff=", 1e-40});
    printfln("mean MPO bond dimension: %f", averageM(H));

    delete[] geom;

    auto psi = run_dmrg(sites, N, sweep_table, sweeps_min, sweeps_max, H, dH2_goal);
    run_analysis(psi);

    return 0;
}
