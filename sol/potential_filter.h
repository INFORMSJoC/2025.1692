// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains all the functions used to manipulate and use
    a dual LP solution in its potential form.
*/

#ifndef POTENTIAL_FILTER_H
#define POTENTIAL_FILTER_H

#include "reduced_obj.h"

namespace HPG {
    struct PotentialFilter {
        // root vertex id
        VertexID root;

        // sum of the duals of the upper bound constraints
        Obj sum_dual;

        // items of the instance
        Items items;

        // table of maximum production of each item in each vertex
        const Table<Cnt>& maximum_production;

        // for each vertex and each item, the corresponding potential
        Table<Obj> potentials;

        // for each vertex its DP value
        Objs dp;

        // for each vertex its RDP value
        Objs rdp;

        // number of cuts
        std::vector<size_t> active_items;

        // create the data structure
        template<typename InstanceT>
        PotentialFilter(const InstanceT& inst, const ReducedObj& obj_func);

        // modify the data structure such that it becomes a copy of the other
        void copy(const PotentialFilter& other);

        // given an item_id, return the reduced profit of the item
        Obj item_obj(int item_id) const;

        // given a hyperarc, return the increased profit of the hyperarc
        Obj hyperarc_obj(const Hyperarc& hyperarc) const;

        // return the improved LP upper bound
        float global_ub() const;

        // return an upper bound on the profit of a solution using the vertex
        float vertex_ub(VertexID r) const;

        // return an upper bound on the profit of a solution using the hyperarc
        float hyperarc_ub(const Hyperarc& hyperarc) const;

        // given a vertex id and an item id, return the potential
        float filter_coeff(VertexID r, int item_id) const;

        // given a vertex id and a label, return the linear part of the profit
        float linear(VertexID r, const std::vector<uint8_t>& cnts) const;

        // given a vertex id, return the constant part of the profit
        float constant(VertexID r) const;

        // given a vertex id and a label, return an upper bound on the profit of a solution using the label
        float label_ub(VertexID r, const std::vector<uint8_t>& cnts) const;
    };
}

namespace HPG {
    template<typename InstanceT>
    PotentialFilter::PotentialFilter(const InstanceT& inst, const ReducedObj& obj_func) 
    :   root(obj_func.root), 
        sum_dual(obj_func.sum_dual), 
        items(obj_func.items), 
        maximum_production(obj_func.maximum_production),
        potentials(inst.nb_vertices(), inst.nb_items(), 0)
    {
        for(int cut_id = 0;cut_id < (int)obj_func.cut_items.size();cut_id++) {
            int item_id = obj_func.cut_items[cut_id];
            for(VertexID r = 0;r < inst.nb_vertices();r++) {
                if(obj_func.activity[r][cut_id / 64] & ((uint64_t)1 << (cut_id % 64))) {
                    potentials.get(r, item_id) += obj_func.cut_mults[cut_id];
                }
            }
        }

        float shifts[items.size()];
        for(int item_id = 0;item_id < (int)items.size();item_id++) {
            shifts[item_id] = items[item_id].obj + potentials.const_get(items[item_id].vrt, item_id);
        }

        for(VertexID r = 0;r < inst.nb_vertices();r++) {
            for(int item_id = 0;item_id < (int)items.size();item_id++) {
                potentials.get(r, item_id) = shifts[item_id] - potentials.get(r, item_id);
            }
        }

        active_items = obj_func.cut_items;
        std::sort(active_items.begin(), active_items.end());
        active_items.erase(std::unique(active_items.begin(), active_items.end()), active_items.end());

        dp = dp_of(inst, *this);
        rdp = rdp_of(inst, *this, dp);
    }

    void PotentialFilter::copy(const PotentialFilter& other) {
        root = other.root;
        sum_dual = other.sum_dual;
        items = other.items;
        potentials = other.potentials;
        dp = other.dp;
        rdp = other.rdp;
        active_items = other.active_items;
    }

    inline Obj PotentialFilter::item_obj(int item_id) const {
        return items[item_id].obj;
    }

    inline Obj PotentialFilter::hyperarc_obj(const Hyperarc& hyperarc) const {
        Obj sum_edge = 0;

        for(size_t item_id : active_items) {
            Obj pr = potentials.const_get(hyperarc.r, item_id);
            Obj pu = potentials.const_get(hyperarc.u, item_id);
            Obj pv = potentials.const_get(hyperarc.v, item_id);

            Obj du = std::max<Obj>(0, pr - pu);
            Obj dv = std::max<Obj>(0, pr - pv);

            Cnt mr = maximum_production.const_get(hyperarc.r, item_id);
            Cnt mu = maximum_production.const_get(hyperarc.u, item_id);
            Cnt mv = maximum_production.const_get(hyperarc.v, item_id);

            if(mu + mv <= mr) {
                sum_edge += mu * du + mv * dv;
            } else if(du > dv) {
                sum_edge += mu * du + std::max(0, mr - mu) * dv;
            } else {
                sum_edge += mv * dv + std::max(0, mr - mv) * du;
            }
        }
        
        if(hyperarc.r == root) return sum_edge + sum_dual;
        return sum_edge;
    }

    inline float PotentialFilter::global_ub() const {
        return dp[root];
    }

    inline float PotentialFilter::vertex_ub(VertexID r) const {
        return dp[r] + rdp[r];
    }

    inline float PotentialFilter::hyperarc_ub(const Hyperarc& hyperarc) const {
        return rdp[hyperarc.r] + dp[hyperarc.u] + dp[hyperarc.v] + hyperarc_obj(hyperarc);
    }

    inline float PotentialFilter::filter_coeff(VertexID r, int item_id) const {
        return potentials.const_get(r, item_id);
    }

    inline float PotentialFilter::linear(VertexID r, const std::vector<uint8_t>& cnts) const {
        Obj sum_objs = 0;
        for(int item_id = 0;item_id < (int)items.size();item_id++) {
            sum_objs += cnts[item_id] * filter_coeff(r, item_id);
        }
        return sum_objs;
    }

    inline float PotentialFilter::constant(VertexID r) const {
        if(r == root) {
            return sum_dual + rdp[r];
        }
        return rdp[r];
    }

    inline float PotentialFilter::label_ub(VertexID r, const std::vector<uint8_t>& cnts) const {
        return constant(r) + linear(r, cnts);
    }
}

#endif