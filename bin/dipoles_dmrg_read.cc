#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "linrot/operators.h"

#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc != 6) {
        printfln("usage: %s <N> <l_max> <sites_path> <H_path> <mps_path>", argv[0]);

        return 1;
    }

    int N = atoi(argv[1]);
    int l_max = atoi(argv[2]);
    std::string sites_path = argv[3];
    std::string H_path = argv[4];
    std::string mps_path = argv[5];

    auto sites = LinearRigidRotor(sites_path, N, {"l_max", l_max});
    auto H = readFromFile<IQMPO>(H_path, sites);
    auto psi = readFromFile<IQMPS>(mps_path, sites);

    printfln("E0 = %.15f", overlap(psi, H, psi));
    run_analysis(psi);

    return 0;
}
