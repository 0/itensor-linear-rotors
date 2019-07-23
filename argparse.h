#ifndef __LINROT_ARGPARSE_H
#define __LINROT_ARGPARSE_H

#include <map>

#include "itensor/all.h"


enum class ArgType {
    Int,
    Real,
    String,
};


class ArgumentTableEntry {
public:
    std::string label;
    ArgType type;
    std::string name;
};


class ArgumentParser {
    std::map<std::string, ArgumentTableEntry> table;

public:
    void add(std::string label, ArgType type) {
        if (table.count(label) > 0) {
            itensor::Error("Argument already present");
        }

        auto name = label;
        name.erase(0, name.find_first_not_of('-'));

        table[label] = ArgumentTableEntry{label, type, name};
    }

    itensor::Args parse(int argc, char* argv[]) {
        std::vector<std::string> arg_strings(argv, argv + argc);
        itensor::Args args;

        for (int idx = 1; idx < argc;) {
            auto current = arg_strings[idx];
            idx++;

            if (table.count(current) == 0) {
                itensor::Error("Unrecognized option \"" + current + "\"");
            }

            auto entry = table[current];

            switch (entry.type) {
                case ArgType::Int:
                    args.add(entry.name, atoi(arg_strings[idx].c_str()));
                    idx++;
                    break;
                case ArgType::Real:
                    args.add(entry.name, atof(arg_strings[idx].c_str()));
                    idx++;
                    break;
                case ArgType::String:
                    args.add(entry.name, arg_strings[idx]);
                    idx++;
                    break;
            }
        }

        return args;
    }
};

#endif // __LINROT_ARGPARSE_H
