#include "itensor/all.h"

#include "linrot/linrigrot.h"
#include "argparse.h"

using namespace itensor;


int main(int argc, char* argv[]) {
    if (argc <= 1) {
        printfln("usage: %s -N <N> --l-max <L> [--lp-sym <S>] [--m-sym <S>]"
                          " --sites-out-path <S>", argv[0]);
        return 1;
    }

    ArgumentParser parser;
    parser.add("-N", ArgType::Int);
    parser.add("--l-max", ArgType::Int);
    parser.add("--lp-sym", ArgType::Bool);
    parser.add("--m-sym", ArgType::Bool);
    parser.add("--sites-out-path", ArgType::String);
    auto args = parser.parse(argc, argv);

    int N = args.getInt("N");
    int l_max = args.getInt("l-max");
    bool lp_sym = args.getBool("lp-sym", true);
    bool m_sym = args.getBool("m-sym", lp_sym);
    auto sites_out_path = args.getString("sites-out-path");

    auto sites = LinearRigidRotor(N, {"l_max", l_max, "lp_sym", lp_sym,
                                      "m_sym", m_sym});

    writeToFile(sites_out_path, sites);

    return 0;
}
