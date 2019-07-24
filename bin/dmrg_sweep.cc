#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "argparse.h"
#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s --sweep-table <T> --num-sweeps <N>"
                          " [--first-sweep <S>] --sites <S> --ham <H>"
                          " [--mps-in <M>] --mps-out <M>", argv[0]);
        return 1;
    }

    ArgumentParser parser;
    parser.add("--sweep-table", ArgType::String);
    parser.add("--num-sweeps", ArgType::Int);
    parser.add("--first-sweep", ArgType::Int, {"required", false});
    parser.add("--sites", ArgType::String);
    parser.add("--ham", ArgType::String);
    parser.add("--mps-in", ArgType::String, {"required", false});
    parser.add("--mps-out", ArgType::String);
    auto args = parser.parse(argc, argv);

    auto sweep_table = InputGroup(args.getString("sweep-table"), "sweeps");
    auto num_sweeps = args.getInt("num-sweeps");
    auto skip_sweeps = args.getInt("first-sweep", 1) - 1;
    auto sites_in_path = args.getString("sites");
    auto H_in_path = args.getString("ham");
    auto mps_in_path = args.getString("mps-in", "");
    auto mps_out_path = args.getString("mps-out");

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
