// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains all the functions that are useful for parsing
    the arguments of the encoder from the command line arguments.
*/

#ifndef ENCODER_ARGUMENTS_H
#define ENCODER_ARGUMENTS_H

#include "stage_grammar.h"
#include "permutation.h"

namespace HPG {
    struct EncoderArguments {
        // dimension of the space
        int dim = 2;

        // for each axis, true if and only if we can use trimming cuts orthogonal to the axis
        std::vector<bool> trim;

        // for each axis, the reduction ratio
        std::vector<int> reduction_ratio;

        // structure of the cutting stages
        StageGrammar stages;
        
        // authorized rotations (ie permutations of the axis)
        std::vector<Permutation> rotations;

        // true if we should round the dimensions up
        // false if we should round the dimensions down
        bool roundup = true;

        // true if the profit of an item is its measure, 
        // false if the original profits are kept
        bool measure = false;

        // true if we should add valid inequalities
        // false otherwise
        bool valid_inequalities = true;

        // path to the input file
        std::string input;
        
        // path to the output file
        std::string output;

        // parse the command line arguments
        EncoderArguments(int argc, char* argv[]);

        // check that the arguments are valid
        void sanitize();
    };
}

namespace HPG {
    EncoderArguments::EncoderArguments(int argc, char* argv[]) {
        std::vector<std::string> args(argv + 1, argv + argc);
        if(args.size() < 2) {
            std::cerr << "missing output file" << std::endl;
            exit(0);
        }

        // iterate over the arguments
        for(size_t i = 0;i < args.size();i++) {
            std::string arg = args[i];

            // split an argument into key and value : 
            // "--dim=3" -> key: "--dim", value: "3"
            auto pos = arg.find('=');
            std::string key = arg.substr(0, pos);
            std::string value = (pos != std::string::npos) ? arg.substr(pos + 1) : "";

            // the equal is optional: for some keys and when value is empty, we parse the next token
            if(value.empty()
            && (key == "-d" || key == "--dim" || key == "-s" || key == "--stages"
            || key == "-t" || key == "--trim" || key == "-r" || key == "--rotations"
            || key == "-h" || key == "--homothety")) {
                if(i + 1 < args.size()) value = args[++i];
                else {
                    std::cerr << "missing value for key " << key << std::endl;
                    exit(0);
                }
            }

            if(key == "-d" || key == "--dim") {
                try {
                    dim = std::stoi(value);
                } catch(const std::exception&) {
                    std::cerr << "the dimension is not an integer" << std::endl;
                    exit(0);
                }
            } else if(key == "-s" || key == "--stages") {
                stages = StageGrammar(value);
            } else if(key == "-t" || key == "--trim") {
                for(char dim_char : value) {
                    trim.push_back(dim_char - '0');
                }
            } else if(key == "-r" || key == "--rotations") {
                Permutation permutation;
                for(char car : value) {
                    if(car == ',') {
                        rotations.push_back(permutation);
                        permutation.clear();
                    }
                    else {
                        permutation.push_back((int)(car - '0'));
                    }
                }
                rotations.push_back(permutation);
            } else if(key == "-h" || key == "--homothety") {
                std::string num;
                for(char car : value) {
                    if(car == ',') {
                        reduction_ratio.push_back(std::stoi(num));
                        num = "";
                    } else {
                        num.push_back(car);
                    }
                }
                reduction_ratio.push_back(std::stoi(num));
            } else if(arg == "-m" || arg == "--measure") {
                measure = true;
            } else if(arg == "-l" || arg == "--low") {
                roundup = false;
            } else if(arg == "-nc" || arg == "--no-cuts") {
                valid_inequalities = false;  
            } else if(i == args.size() - 2) {
                input = arg;
            } else if(i == args.size() - 1) {
                output = arg;
            } else {
                std::cerr << "unknown key " << key << std::endl;
                exit(0);
            }
        }

        // initialize with default values

        // homothety by default
        if(reduction_ratio.empty()) {
            reduction_ratio.resize(dim, 1);
        }

        // any-stage by default
        if(stages.states.empty()) {
            stages = StageGrammar();
            stages.states.push_back({std::vector<int>(dim, 0)});
            std::iota(stages.states[0].axes.begin(), stages.states[0].axes.end(), 0);
            stages.start = stages.end = 0;
        }

        // allow trimming by default
        if(trim.empty()) {
            trim.resize(dim, true);
        }
        
        // all permutations by default
        if(rotations.empty()) {
            rotations = all_permutations(dim);
        }

        // check that the arguments are valid
        sanitize();
    }

    void EncoderArguments::sanitize() {
        // dim check
        if(dim <= 0 || dim > 10) {
            std::cerr << "invalid --dim value" << std::endl;
            exit(0);
        }

        // stages check
        for(const auto& state : stages.states) {
            for(const auto& axis : state.axes) {
                if(axis < 0 || axis >= dim) {
                    std::cout << "invalid --stages value" << std::endl;
                    exit(0);
                }
            }
        }

        // trim check
        if(trim.size() != dim) {
            std::cout << "invalid --trim value" << std::endl;
            exit(0);
        }

        // rotations check
        for(const auto& permutation : rotations) {
            bool valid = permutation.size() == dim;

            std::vector<bool> included(dim, false);
            for(int axis : permutation) {
                valid = valid && axis >= 0 && axis < dim && !included[axis];
                included[axis] = true;
            }

            if(!valid) {
                std::cout << "invalid --permutation value" << std::endl;
                exit(0);
            }
        }
    }
}

#endif