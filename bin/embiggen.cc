#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "argparse.h"
#include "dmrg.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s --sites1 <S> --mps1-in <M> --sites2 <S>"
                          " --mps2-out <M>", argv[0]);
        return 1;
    }

    ArgumentParser parser;
    parser.add("--sites1", ArgType::String);
    parser.add("--mps1-in", ArgType::String);
    parser.add("--sites2", ArgType::String);
    parser.add("--mps2-out", ArgType::String);
    auto args = parser.parse(argc, argv);

    auto sites1_in_path = args.getString("sites1");
    auto mps1_in_path = args.getString("mps1-in");
    auto sites2_in_path = args.getString("sites2");
    auto mps2_out_path = args.getString("mps2-out");

    auto sites1 = readFromFile<LinearRigidRotor>(sites1_in_path);
    auto mps1 = readFromFile<MPS>(mps1_in_path, sites1);
    auto sites2 = readFromFile<LinearRigidRotor>(sites2_in_path);

    if (length(sites1) != length(sites2)) {
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
