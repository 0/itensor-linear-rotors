#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "argparse.h"
#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s --sweep-table <T> [--mps-cutoff <C>]"
                          " --num-sweeps <N> [--first-sweep <S>] --sites <S>"
                          " --ham <H> [--ortho-mps <M>] [--mps-in <M>]"
                          " --mps-out <M>", argv[0]);
        return 1;
    }

    ArgumentParser parser;
    parser.add("--sweep-table", ArgType::String);
    parser.add("--mps-cutoff", ArgType::Real, {"required", false});
    parser.add("--num-sweeps", ArgType::Int);
    parser.add("--first-sweep", ArgType::Int, {"required", false});
    parser.add("--sites", ArgType::String);
    parser.add("--ham", ArgType::String);
    parser.add("--ortho-mps", ArgType::String, {"required", false});
    parser.add("--mps-in", ArgType::String, {"required", false});
    parser.add("--mps-out", ArgType::String);
    auto args = parser.parse(argc, argv);

    auto sweep_table = InputGroup(args.getString("sweep-table"), "sweeps");
    auto mps_cutoff = args.getReal("mps-cutoff", -1.0);
    auto num_sweeps = args.getInt("num-sweeps");
    auto skip_sweeps = args.getInt("first-sweep", 1) - 1;
    auto sites_in_path = args.getString("sites");
    auto H_in_path = args.getString("ham");
    auto ortho_mps_in_path = args.getString("ortho-mps", "");
    auto mps_in_path = args.getString("mps-in", "");
    auto mps_out_path = args.getString("mps-out");

    auto sites = readFromFile<LinearRigidRotor>(sites_in_path);
    auto N = length(sites);
    auto H = readFromFile<MPO>(H_in_path, sites);

    std::vector<MPS> ortho_wfs;
    if (!ortho_mps_in_path.empty()) {
        ortho_wfs.push_back(readFromFile<MPS>(ortho_mps_in_path, sites));
    }

    MPS psi;
    if (mps_in_path.empty()) {
        auto state = InitState(sites);
        for (auto i : range1(N)) {
            state.set(i, "l0m0");
        }
        psi = MPS(state);
    } else {
        psi = readFromFile<MPS>(mps_in_path, sites);
    }

    cpu_time time;
    dmrg_sweep(psi, H, sweep_table, num_sweeps, skip_sweeps, ortho_wfs, mps_cutoff);
    printfln("cputime = %.15e", time.sincemark().time);
    printfln("walltime = %.15e", time.sincemark().wall);

    writeToFile(mps_out_path, psi);

    return 0;
}
