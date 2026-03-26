// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains the column and row generation algorithm to compute
    linear programming relaxations.
*/

#ifndef LP_H
#define LP_H

#include "dp.h"
#include "model.h"
#include "potential_filter.h"

namespace HPG {
    // threshold for which a cut is considered violated
    const float CUT_TOLERANCE = 1e-9;

    // value for the convex combination between the current solution and the smoothed solution
    const float STABILISATION = 0.2;

    template<typename InstanceT>
    struct LpThread {
        const InstanceT& inst;
        Table<Cnt> maximum_production;
        Model model;
        LpSolution smoothed_solution;
        bool must_stop;

        // compute the optimal solution without cuts
        LpThread(const InstanceT& inst);

        // reoptimize the model using column (hyperarc) generation
        PotentialFilter hyperarc_generation();

        // add item flow cuts to the model and reoptimize
        PotentialFilter cut_generation();
    };
}

namespace HPG {
    template<typename InstanceT>
    LpThread<InstanceT>::LpThread(const InstanceT& _inst)
    :   inst(_inst),
        maximum_production(maximum_production_of(inst)),
        model(inst, maximum_production),
        must_stop(false)
    {
        smoothed_solution.item_flow.resize(inst.nb_items(), 0);
    }

    template<typename InstanceT>
    PotentialFilter LpThread<InstanceT>::hyperarc_generation() {
        logger.nb_lp_iters = 0;

        while(!must_stop) {
            logger.nb_lp_iters += 1;
            model.optimize(true);

            auto reduced_obj = ReducedObj(inst, model);
            auto dp = dp_of(inst, reduced_obj);
            auto sol = dp_solution_of(inst, reduced_obj, dp);

            auto hyperarcs_before = model.nb_hyperarcs();
            sol.enum_up([&] (VertexID r) {
                sol.enum_hyperarcs(r, [&] (const Hyperarc& edge) {
                    model.new_hyperarc(edge);
                });
            });
            auto hyperarcs_after = model.nb_hyperarcs();

            if(hyperarcs_after == hyperarcs_before) break;
        }

        if(!must_stop) {
            auto current_solution = model.get_solution();
            logger.nb_hyperarcs = model.nb_hyperarcs();
            logger.active_hyperarcs = current_solution.hyperarc_flow.size();
            logger.nb_cuts = model.cut_sets.size();
        }

        return PotentialFilter(inst, ReducedObj(inst, model));
    }

    template<typename InstanceT>
    PotentialFilter LpThread<InstanceT>::cut_generation() {
        auto current_solution = model.get_solution();
        smoothed_solution = (1 - STABILISATION) * smoothed_solution + STABILISATION * current_solution;

        for(int item_id = 0;item_id < inst.nb_items();item_id++) {
            if(smoothed_solution.item_flow[item_id] > CUT_TOLERANCE) {
                auto cut = cut_generator(
                    inst.nb_vertices(), 
                    item_id, 
                    inst.items[item_id].vrt, 
                    inst.root, 
                    smoothed_solution, 
                    maximum_production
                );
                if(cut.value + CUT_TOLERANCE > smoothed_solution.item_flow[item_id]) continue;
                model.new_cut(item_id, cut.vertex_side);
            }
        }
        model.optimize(false);

        return hyperarc_generation();
    }
}

#endif