#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "linrot/operators.h"
#include "argparse.h"

#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s --sites-path <S> --H-path <H>"
                          " --mps-path <M>", argv[0]);
        return 0;
    }

    ArgumentParser parser;
    parser.add("--sites-path", ArgType::String);
    parser.add("--H-path", ArgType::String);
    parser.add("--mps-path", ArgType::String);
    auto args = parser.parse(argc, argv);

    auto sites_path = args.getString("sites-path");
    auto H_path = args.getString("H-path");
    auto mps_path = args.getString("mps-path");

    auto sites = readFromFile<LinearRigidRotor>(sites_path);
    auto H = readFromFile<IQMPO>(H_path, sites);
    auto psi = readFromFile<IQMPS>(mps_path, sites);

    printfln("E0 = %.15f", overlap(psi, H, psi));
    run_analysis(psi);

    return 0;
}
