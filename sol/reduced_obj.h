// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains all the functions used to manipulate and use
    a dual LP solution in its raw form.
*/

#ifndef REDUCED_OBJ_H
#define REDUCED_OBJ_H

#include "model.h"

namespace HPG {
    using Bitset = std::vector<uint64_t>;

    struct ReducedObj {
        // root vertex id
        VertexID root;

        // sum of the duals of the upper bound constraints
        Obj sum_dual;

        // items of the instance
        Items items;

        // table of maximum production of each item in each vertex
        const Table<Cnt>& maximum_production;

        // for each cut, its item
        std::vector<size_t> cut_items;

        // for each cut, its dual value
        std::vector<float> cut_mults;

        // for each vertex, a bitset representing the cuts for which the vertex is in the cut's set
        std::vector<Bitset> activity;

        // create the data structure
        template<typename InstanceT>
        ReducedObj(const InstanceT& inst, const Model& model);

        // given an item_id, return the reduced profit of the item
        Obj item_obj(int item_id) const;

        // given a hyperarc, return the increased profit of the hyperarc
        Obj hyperarc_obj(const Hyperarc& hyperarc) const;
    };
}

namespace HPG {
    template<typename InstanceT>
    ReducedObj::ReducedObj(const InstanceT& inst, const Model& model) 
    :   root(inst.root), 
        sum_dual(0), 
        items(inst.items), 
        maximum_production(model.maximum_production) 
    {
        for(int item_id = 0;item_id < (int)items.size();item_id++) {
            items[item_id].obj += model.get_item_dual(item_id);
            sum_dual += std::max(0., items[item_id].ub * -model.get_item_dual(item_id));
        }

        std::vector<size_t> active_cuts;
        for(int cut_id = 0;cut_id < (int)model.cut_sets.size();cut_id++) {
            double dual = -model.get_cut_dual(cut_id);
            items[model.cut_items[cut_id]].obj -= dual;

            if(dual > 0) {
                cut_items.push_back(model.cut_items[cut_id]);
                cut_mults.push_back(dual);
                active_cuts.push_back(cut_id);
            }
        }

        activity.resize(inst.nb_vertices(), Bitset((active_cuts.size() + 63) / 64, 0));
        for(VertexID r = 0;r < inst.nb_vertices();r++) {
            for(int active_id = 0;active_id < active_cuts.size();active_id++) {
                int cut_id = active_cuts[active_id];
                if(model.cut_sets[cut_id][r]) {
                    activity[r][active_id / 64] |= ((uint64_t)1 << (active_id % 64));
                }
            }
        }
    }

    inline Obj ReducedObj::item_obj(int item_id) const {
        return items[item_id].obj;
    }

    inline Obj ReducedObj::hyperarc_obj(const Hyperarc& hyperarc) const {
        Obj sum_edge = 0;

        for(int chunk_id = 0;chunk_id < activity[hyperarc.r].size();chunk_id++) {
            uint64_t bits = ~activity[hyperarc.r][chunk_id] & (activity[hyperarc.u][chunk_id] | activity[hyperarc.v][chunk_id]);

            while(bits != 0) {
                int bit_id = std::countr_zero(bits);
                int cut_id = 64 * chunk_id + bit_id;
                sum_edge += cut_coeff(cut_items[cut_id], hyperarc, 
                    activity[hyperarc.u][chunk_id] & ((uint64_t)1 << bit_id), 
                    activity[hyperarc.v][chunk_id] & ((uint64_t)1 << bit_id), 
                    maximum_production
                ) * cut_mults[cut_id];
                bits &= bits - 1;
            }
        }

        if(hyperarc.r == root) return sum_edge + sum_dual;
        return sum_edge;
    }
}

#endif