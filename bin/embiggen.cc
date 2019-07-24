#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "argparse.h"
#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s --sites1-in-path <S> --mps1-in-path <M>"
                          " --sites2-in-path <S>"
                          " --mps2-out-path <M>", argv[0]);
        return 1;
    }

    ArgumentParser parser;
    parser.add("--sites1-in-path", ArgType::String);
    parser.add("--mps1-in-path", ArgType::String);
    parser.add("--sites2-in-path", ArgType::String);
    parser.add("--mps2-out-path", ArgType::String);
    auto args = parser.parse(argc, argv);

    auto sites1_in_path = args.getString("sites1-in-path");
    auto mps1_in_path = args.getString("mps1-in-path");
    auto sites2_in_path = args.getString("sites2-in-path");
    auto mps2_out_path = args.getString("mps2-out-path");

    auto sites1 = readFromFile<LinearRigidRotor>(sites1_in_path);
    auto mps1 = readFromFile<IQMPS>(mps1_in_path, sites1);
    auto sites2 = readFromFile<LinearRigidRotor>(sites2_in_path);

    if (sites1.N() != sites2.N()) {
        println("Number of sites must match");
        return 1;
    }
    if (sites1.l_max() >= sites2.l_max()) {
        println("Local basis must be larger");
        return 1;
    }

    auto mps2 = embiggen(sites1, mps1, sites2);

    writeToFile(mps2_out_path, mps2);

    return 0;
}
