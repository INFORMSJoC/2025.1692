// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains all the functions related to directed graphs.
*/

#ifndef DIGRAPH_H
#define DIGRAPH_H

#include "../common/types.h"

namespace HPG {
    struct DigraphCut {
        // the value of the cut
        double value;

        // for each vertex, tell if the vertex and the source have the same side 
        std::vector<bool> vertex_side;
    };

    struct DigraphArc {
        // id of the reverse arc
        int reverse_arc_id;

        // capacity of the arc
        double capacity;
    };

    struct DigraphNeighbor {
        // id of the neighbor
        VertexID vertex;

        // id of the arc going to the vertex
        int arc_id;
    };

    // threshold for considering a capacity to be zero
    const double DINIC_PRECISION = 1e-9;

    struct Digraph {
        // contain every arc in the digraph
        std::vector<DigraphArc> arcs;

        // for each vertex, contain all of its neighbors
        std::vector<std::vector<DigraphNeighbor>> neighbors;

        // create a new vertex and return it
        VertexID new_vertex();

        // create a new arc form u to v with a given capacity
        void new_arc(VertexID u, VertexID v, double capacity);

        // compute the minimum cut between the vertices source and sink using Dinic's algorithm
        DigraphCut min_cut(VertexID source, VertexID sink);
    };
}

namespace HPG {
    VertexID Digraph::new_vertex() {
        VertexID id = neighbors.size();
        neighbors.push_back({});
        return id;
    }

    void Digraph::new_arc(VertexID u, VertexID v, double capacity) {
        int arc_id = arcs.size();
        int reverse_arc_id = arc_id + 1;
        arcs.push_back({reverse_arc_id, capacity});
        arcs.push_back({arc_id, 0.});
        neighbors[u].push_back({v, arc_id});
        neighbors[v].push_back({u, reverse_arc_id});
    }

    DigraphCut Digraph::min_cut(VertexID source, VertexID sink) {
        const int MAXIMUM_DISTANCE = 1 + neighbors.size();

        double cut_value = 0;

        // for each arc, the flow that traverse it
        // the value is negative for the reverse arc
        std::vector<double> arc_flow(arcs.size(), 0);

        // for each vertex, its distance from the source
        std::vector<int> distances(neighbors.size());

        // compute the residual capacity of an arc
        auto residual_capacity = [&] (int arc_id) {
            return arcs[arc_id].capacity - arc_flow[arc_id];
        };

        // recompute the distances from the source in the directed graph with only the residual arcs using BFS
        auto recompute_distances = [&] () {
            fill(distances.begin(), distances.end(), MAXIMUM_DISTANCE);
            distances[source] = 0;

            std::deque<int> active_vertices = { source };
            while(!active_vertices.empty()) {
                int u = active_vertices.front();
                active_vertices.pop_front();

                for(auto neighbor : neighbors[u]) {
                    if(distances[neighbor.vertex] == MAXIMUM_DISTANCE 
                    && residual_capacity(neighbor.arc_id) > DINIC_PRECISION) {
                        distances[neighbor.vertex] = distances[u] + 1;
                        active_vertices.push_back(neighbor.vertex);
                    }
                }
            }
        };

        // for each vertex, the first position of a candidate arc in its neighbors  
        std::vector<int> current_arc_pos(neighbors.size());

        // try to push a at most available_flow from u to the sink, and return the pushed flow
        auto push_flow = [&] (auto self, VertexID u, double available_flow) -> double {
            if(u == sink) return available_flow;

            while(current_arc_pos[u] != neighbors[u].size()) {
                DigraphNeighbor& neighbor = neighbors[u][current_arc_pos[u]];
                if(distances[neighbor.vertex] == distances[u] + 1
                && residual_capacity(neighbor.arc_id) > DINIC_PRECISION) {
                    double pushed_flow = self(self, 
                        neighbor.vertex, 
                        std::min(available_flow, residual_capacity(neighbor.arc_id))
                    );

                    if(pushed_flow > DINIC_PRECISION) {
                        arc_flow[neighbor.arc_id] += pushed_flow;
                        arc_flow[arcs[neighbor.arc_id].reverse_arc_id] -= pushed_flow;
                        return pushed_flow;
                    }
                }
                current_arc_pos[u]++;
            }

            return 0;
        };

        // main loop of Dinic's algorithm
        do {
            recompute_distances();
            std::fill(current_arc_pos.begin(), current_arc_pos.end(), 0);

            // compute a blocking flow with shortest paths
            while(true) {
                double pushed_flow = push_flow(push_flow, source, INFINITY);
                if(pushed_flow < DINIC_PRECISION) break;
                else cut_value += pushed_flow;
            }
        } while(distances[sink] != MAXIMUM_DISTANCE);


        // recursively explore the digraph to determine the side of each vertex
        DigraphCut cut = { cut_value, std::vector<bool>(neighbors.size(), false) };
        auto mark_source_side_rec = [&] (auto self, int u) -> void {
            if(!cut.vertex_side[u]) {
                cut.vertex_side[u] = true;

                for(auto neighbor : neighbors[u]) {
                    if(residual_capacity(neighbor.arc_id) > DINIC_PRECISION) {
                        self(self, neighbor.vertex);
                    }
                }
            }
        };
        mark_source_side_rec(mark_source_side_rec, source);

        return cut;
    }
}

#endif