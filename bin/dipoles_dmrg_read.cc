#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "linrot/operators.h"

#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc != 4) {
        printfln("usage: %s <sites_path> <H_path> <mps_path>", argv[0]);

        return 1;
    }

    std::string sites_path = argv[1];
    std::string H_path = argv[2];
    std::string mps_path = argv[3];

    auto sites = readFromFile<LinearRigidRotor>(sites_path);
    auto H = readFromFile<IQMPO>(H_path, sites);
    auto psi = readFromFile<IQMPS>(mps_path, sites);

    printfln("E0 = %.15f", overlap(psi, H, psi));
    run_analysis(psi);

    return 0;
}
