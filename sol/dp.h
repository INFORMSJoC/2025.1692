// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains the dynamic programming algorithms on hypergraphs : 
    all these algorithms ignore the upper bounds of the items.
*/

#ifndef DP_H
#define DP_H

#include "table.h"

namespace HPG {
    // in what follows, an objective function (ObjFuncT) is an object with two methods :
    // hyperarc_obj : given a hyperarc, return its profit
    // item_obj : given an item id, return its profit

    // given a hypergraph and an objective function, return an array :
    // for each vertex, the maximum profit achievable by a decomposition of the vertex.
    template<typename InstanceT, typename ObjFuncT>
    Objs dp_of(const InstanceT& inst, const ObjFuncT& obj_func);

    // given a hypergraph, an objective function and the result of the dp_of function, return an array :
    // for each vertex, the maximum profit achievable by a decomposition of the root giving the vertex.
    template<typename InstanceT, typename ObjFuncT>
    Objs rdp_of(const InstanceT& inst, const ObjFuncT& obj_func, const Objs& dp);

    // given a hypergraph, an objective function and the result of the dp_of function, return a hypergraph
    // which represents the decomposition of the root of maximum profit.
    template<typename InstanceT, typename ObjFuncT>
    HPG dp_solution_of(const InstanceT& inst, const ObjFuncT& obj_func, const Objs& dp);

    // given a hypergraph, return a table : for each vertex and each item
    // the maximum number of times the item can be produced by a decomposition of the vertex
    template<typename InstanceT>
    Table<Cnt> maximum_production_of(const InstanceT& inst);
};

namespace HPG {
    template<typename InstanceT, typename ObjFuncT>
    Objs dp_of(const InstanceT& inst, const ObjFuncT& obj_func) {
        Objs dp(inst.nb_vertices(), 0);
        inst.enum_up([&] (VertexID r) {
            for(int item_id = 0;item_id < inst.nb_items();item_id++) {
                if(inst.items[item_id].vrt == r) {
                    dp[r] = std::max(dp[r], obj_func.item_obj(item_id));
                }
            }

            inst.enum_hyperarcs(r, [&](const Hyperarc& h) {
                dp[r] = std::max(dp[r], dp[h.u] + dp[h.v] + obj_func.hyperarc_obj(h));
            });
        });
        return dp;
    }

    template<typename InstanceT, typename ObjFuncT>
    Objs rdp_of(const InstanceT& inst, const ObjFuncT& obj_func, const Objs& dp) {
        Objs rdp(inst.nb_vertices(), -INFINITY);
        rdp[inst.root] = 0;

        inst.enum_down([&] (VertexID r) {
            inst.enum_hyperarcs(r, [&] (const Hyperarc& h) {
                Obj obj = obj_func.hyperarc_obj(h);
                rdp[h.u] = std::max(rdp[h.u], rdp[h.r] + dp[h.v] + obj);
                rdp[h.v] = std::max(rdp[h.v], rdp[h.r] + dp[h.u] + obj);
            });
        });

        return rdp;
    }

    template<typename InstanceT, typename ObjFuncT>
    HPG dp_solution_of(const InstanceT& inst, const ObjFuncT& obj_func, const Objs& dp) {
        HPG sol = remove_hyperarcs(inst);

        auto create_solution_rec = [&](auto self, int r) -> void {
            if(dp[r] == 0) return;

            for(int item_id = 0;item_id < inst.nb_items();item_id++) {
                if(inst.items[item_id].vrt == r && obj_func.item_obj(item_id) == dp[r]) {
                    return;
                }
            }

            Hyperarc hyperarc;
            inst.enum_hyperarcs(r, [&](const Hyperarc& h) {
                if(dp[h.u] + dp[h.v] + obj_func.hyperarc_obj(h) == dp[h.r]) {
                    hyperarc = h;
                }
            });

            sol.new_hyperarc(hyperarc);
            self(self, hyperarc.u);
            self(self, hyperarc.v);
        };

        create_solution_rec(create_solution_rec, inst.root);
        return sol;
    }

    template<typename InstanceT>
    Table<Cnt> maximum_production_of(const InstanceT& inst) {
        Table<Cnt> table(inst.nb_vertices(), inst.nb_items(), 0);

        inst.enum_up([&] (VertexID r) {
            inst.enum_hyperarcs(r, [&](const Hyperarc& h) {
                for(int item_id = 0;item_id < inst.nb_items();item_id++) {
                    table.get(h.r, item_id) = std::max(
                        table.get(h.r, item_id),
                        table.get(h.u, item_id) + table.get(h.v, item_id)
                    );
                }
            });

            for(int item_id = 0;item_id < inst.nb_items();item_id++) {
                if(inst.items[item_id].vrt == r) {
                    table.get(r, item_id) = std::max(table.get(r, item_id), 1);
                }
                table.get(r, item_id) = std::min(table.get(r, item_id), inst.items[item_id].ub);
            }
        });

        return table;
    }
};

#endif