#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "argparse.h"
#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s --sites <S> --mps <M>", argv[0]);
        return 1;
    }

    ArgumentParser parser;
    parser.add("--sites", ArgType::String);
    parser.add("--mps", ArgType::String);
    auto args = parser.parse(argc, argv);

    auto sites_in_path = args.getString("sites");
    auto mps_in_path = args.getString("mps");

    auto sites = readFromFile<LinearRigidRotor>(sites_in_path);
    auto psi = readFromFile<MPS>(mps_in_path, sites);

    run_analysis(sites, psi);

    return 0;
}
