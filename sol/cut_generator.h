// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains the separation algorithm for the item flow cuts.
*/

#ifndef CUT_GENERATOR_H
#define CUT_GENERATOR_H

#include "digraph.h"
#include "lp_solution.h"
#include "table.h"
#include "logger.h"

namespace HPG {
    // given an item id, an hyperarc (r, u, v), the sides of u and v in a cut of the item, and the maximum production table
    // return the coefficient of the variable of the hyperarc in the cut's inequality.
    Cnt cut_coeff(int item_id, const Hyperarc& hyperarc, bool side_u, bool side_v, const Table<Cnt>& maximum_production);

    // generate a cut from a fractional solution
    DigraphCut cut_generator(int nb_vertices, int item_id, VertexID source, VertexID sink, const LpSolution& sol, const Table<Cnt>& maximum_production);
}

namespace HPG {
    inline Cnt cut_coeff(int item_id, const Hyperarc& hyperarc, bool side_u, bool side_v, const Table<Cnt>& maximum_production) {
        Cnt max_u = maximum_production.const_get(hyperarc.u, item_id);
        Cnt max_v = maximum_production.const_get(hyperarc.v, item_id);
        Cnt max_e = std::min(max_u + max_v, maximum_production.const_get(hyperarc.r, item_id));
        
        if(side_u && side_v) return max_e;
        else if(side_u) return max_u;
        else if(side_v) return max_v;
        else return 0;
    }

    DigraphCut cut_generator(int nb_vertices, int item_id, VertexID source, VertexID sink, const LpSolution& sol, const Table<Cnt>& maximum_production) {
        Digraph g;
        std::map<VertexID, VertexID> compression;

        auto get_vertex = [&] (int vertex) {
            auto it_vertex = compression.find(vertex);
            if(it_vertex != compression.end()) return it_vertex->second;
            VertexID id = g.new_vertex();
            return compression[vertex] = id;
        };

        // create the digraph associated with the fractional solution
        for(auto it : sol.hyperarc_flow) {
            const Hyperarc& h = it.first;
            double x = it.second;
            VertexID hub = g.new_vertex();

            g.new_arc(get_vertex(h.u), hub, x * cut_coeff(item_id, h, true, false, maximum_production));
            g.new_arc(get_vertex(h.v), hub, x * cut_coeff(item_id, h, false, true, maximum_production));
            g.new_arc(hub, get_vertex(h.r), x * cut_coeff(item_id, h, true, true, maximum_production));
        }

        // compute the minimum cut in the digraph
        auto cut = g.min_cut(get_vertex(source), get_vertex(sink));
        
        // compute the side of each vertex in the hypergraph
        std::vector<bool> hpg_side(nb_vertices, false);
        for(auto it : compression) {
            hpg_side[it.first] = cut.vertex_side[it.second];
        }

        return { cut.value, hpg_side };
    }
}

#endif