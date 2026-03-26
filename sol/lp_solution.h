// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains the implementation of the operators 
    for manipulating a LP solution. 
*/

#ifndef LP_SOLUTION_H
#define LP_SOLUTION_H

#include "../common/hpg.h"

namespace HPG {
    struct LpSolution {
        // for each item, the number of times it is produced
        std::vector<double> item_flow;

        // for each hyperarc, the numer of times it is used
        std::map<Hyperarc, double> hyperarc_flow;
    };

    // multiplication of a scalar by a fractional solution 
    LpSolution operator* (const double& scalar, const LpSolution& solution);

    // addition of two fractional solutions
    LpSolution operator+ (const LpSolution& solution_a, const LpSolution& solution_b);
}

namespace HPG {
    LpSolution operator* (const float& scalar, const LpSolution& solution) {
        LpSolution result;
        result.item_flow.resize(solution.item_flow.size());
        for(int item_id = 0;item_id < (int)solution.item_flow.size();item_id++) {
            result.item_flow[item_id] = scalar * solution.item_flow[item_id];
        }

        if(scalar != 0) {
            for(auto hyperarc_iterator : solution.hyperarc_flow) {
                result.hyperarc_flow[hyperarc_iterator.first] = scalar * hyperarc_iterator.second;
            }
        }

        return result;
    }

    LpSolution operator+ (const LpSolution& solution_a, const LpSolution& solution_b) {
        LpSolution result;
        result.item_flow.resize(solution_a.item_flow.size());
        for(int item_id = 0;item_id < (int)solution_a.item_flow.size();item_id++) {
            result.item_flow[item_id] = solution_a.item_flow[item_id] + solution_b.item_flow[item_id];
        }

        for(auto hyperarc_iterator : solution_a.hyperarc_flow) {
            result.hyperarc_flow[hyperarc_iterator.first] += hyperarc_iterator.second;
        }

        for(auto hyperarc_iterator : solution_b.hyperarc_flow) {
            result.hyperarc_flow[hyperarc_iterator.first] += hyperarc_iterator.second;
        }

        return result;
    }
}

#endif