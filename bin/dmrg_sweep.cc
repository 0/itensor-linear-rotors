#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "argparse.h"
#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s --sweep-table <T> --num-sweeps <N>"
                          " [--skip-sweeps <S>] --sites-in-path <S>"
                          " --H-in-path <H> [--mps-in-path <M>]"
                          " --mps-out-path <M>", argv[0]);
        return 1;
    }

    ArgumentParser parser;
    parser.add("--sweep-table", ArgType::String);
    parser.add("--num-sweeps", ArgType::Int);
    parser.add("--skip-sweeps", ArgType::Int, {"required", false});
    parser.add("--sites-in-path", ArgType::String);
    parser.add("--H-in-path", ArgType::String);
    parser.add("--mps-in-path", ArgType::String, {"required", false});
    parser.add("--mps-out-path", ArgType::String);
    auto args = parser.parse(argc, argv);

    auto sweep_table = InputGroup(args.getString("sweep-table"), "sweeps");
    auto num_sweeps = args.getInt("num-sweeps");
    auto skip_sweeps = args.getInt("skip-sweeps", 0);
    auto sites_in_path = args.getString("sites-in-path");
    auto H_in_path = args.getString("H-in-path");
    auto mps_in_path = args.getString("mps-in-path", "");
    auto mps_out_path = args.getString("mps-out-path");

    auto sites = readFromFile<LinearRigidRotor>(sites_in_path);
    auto N = sites.N();
    auto H = readFromFile<IQMPO>(H_in_path, sites);

    IQMPS psi;
    if (mps_in_path.empty()) {
        auto state = InitState(sites);
        for (auto i : range1(N)) {
            state.set(i, "l0m0");
        }
        psi = IQMPS(state);
    } else {
        psi = readFromFile<IQMPS>(mps_in_path, sites);
    }

    dmrg_sweep(psi, H, sweep_table, num_sweeps, skip_sweeps);

    writeToFile(mps_out_path, psi);

    return 0;
}
