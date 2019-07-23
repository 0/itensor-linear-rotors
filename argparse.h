#ifndef __LINROT_ARGPARSE_H
#define __LINROT_ARGPARSE_H

#include <map>

#include "itensor/all.h"


enum class ArgType {
    Flag,
    Int,
    Bool,
    Real,
    String,
};


class ArgumentTableEntry {
public:
    std::string label;
    ArgType type;
    std::string name;
};


template<class T>
std::set<T> set_diff(std::set<T> a, std::set<T> b) {
    std::set<T> amb;
    std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(amb, amb.end()));
    return amb;
}


class ArgumentParser {
    std::map<std::string, ArgumentTableEntry> table;
    std::set<std::string> required;

public:
    void add(std::string label, ArgType type, itensor::Args const& args = itensor::Args::global()) {
        if (table.count(label) > 0) {
            itensor::Error("Argument already present");
        }

        auto name = label;
        name.erase(0, name.find_first_not_of('-'));

        table[label] = ArgumentTableEntry{label, type, name};

        if (args.getBool("required", true)) {
            required.insert(label);
        }
    }

    itensor::Args parse(int argc, char* argv[]) {
        std::vector<std::string> arg_strings(argv, argv + argc);
        std::set<std::string> seen;
        itensor::Args args;

        for (int idx = 1; idx < argc;) {
            auto current = arg_strings[idx];
            idx++;

            if (table.count(current) == 0) {
                itensor::Error("Unrecognized option \"" + current + "\"");
            }

            auto entry = table[current];
            seen.insert(current);

            switch (entry.type) {
                case ArgType::Flag:
                    args.add(entry.name, true);
                    break;
                case ArgType::Int:
                    args.add(entry.name, atoi(arg_strings[idx].c_str()));
                    idx++;
                    break;
                case ArgType::Bool:
                    args.add(entry.name, static_cast<bool>(atoi(arg_strings[idx].c_str())));
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

        auto missing = set_diff(required, seen);
        if (missing.size() > 0) {
            std::ostringstream oss;
            oss << "Required options not given:";
            for (auto it = missing.begin(); it != missing.end(); it++) {
                oss << " " << *it;
            }
            itensor::Error(oss.str());
        }

        return args;
    }
};

#endif // __LINROT_ARGPARSE_H
