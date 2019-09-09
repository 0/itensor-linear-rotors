#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "argparse.h"
#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s --num-samples <N> --sites-in-path <S>"
                          " --mps-in-path <M>", argv[0]);
        return 1;
    }

    ArgumentParser parser;
    parser.add("--num-samples", ArgType::Int);
    parser.add("--sites-in-path", ArgType::String);
    parser.add("--mps-in-path", ArgType::String);
    auto args = parser.parse(argc, argv);

    int num_samples = args.getInt("num-samples");
    auto sites_in_path = args.getString("sites-in-path");
    auto mps_in_path = args.getString("mps-in-path");

    auto sites = readFromFile<LinearRigidRotor>(sites_in_path);
    auto psi = readFromFile<IQMPS>(mps_in_path, sites);

    run_sampling(sites.l_max(), psi, num_samples);

    return 0;
}
