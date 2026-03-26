// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains the basic algorithms on hypergraphs.
*/

#ifndef HPG_ALGORITHMS_H
#define HPG_ALGORITHMS_H

#include "hpg.h"

namespace HPG {
    // return the total number of hyperarcs of the hypergraph
    template<typename InstanceT>
    size_t nb_hyperarcs(const InstanceT& inst);

    // save the hypergraph to a hpg file
    template<typename InstanceT>
    void save_hpg(const InstanceT& inst, const std::string& path);

    // return for each item, the number of times its vertex is the source of an hyperarc
    template<typename InstanceT>
    Cnts usage_of_items(const InstanceT& inst);

    // return the total profit of the given hypergraph, interpreted as a solution
    template<typename InstanceT>
    Obj solution_value(const InstanceT& inst);

    // compute the topological ordering of the vertices
    template<typename InstanceT>
    std::vector<VertexID> compute_order(const InstanceT& inst);
}

namespace HPG {
    template<typename InstanceT>
    size_t nb_hyperarcs(const InstanceT& inst) {
        size_t nb = 0;
        inst.enum_up([&] (VertexID r) {
            inst.enum_hyperarcs(r, [&] (const Hyperarc& _) {
                nb += 1;
            });
        });
        return nb;
    }

    template<typename InstanceT>
    void save_hpg(const InstanceT& inst, const std::string& path) {
        std::ofstream fout;
        fout.open(path);

        if(!fout.is_open()) {
            std::cout << "unable to open " << path << std::endl;
            exit(0);
        }

        fout << "<symbols>\n";
        for(const auto& name : inst.vertex_names) {
            fout << name << "\n";
        }
        fout << "<edges>\n";
        for(VertexID r = 0;r < inst.nb_vertices();r++) {
            fout << r << "\n";
            inst.enum_hyperarcs(r, [&] (Hyperarc hyperarc) {
                fout << hyperarc.u << " " << hyperarc.v << "\n";
            });
        };
        fout << "<root>\n";
        fout << inst.root << "\n";
        fout << "<items>\n";
        for(const auto& item : inst.items) {
            fout << item.vrt << " " << item.obj << " " << item.ub << "\n";
        }

        fout.close();
    }

    template<typename InstanceT>
    Cnts usage_of_items(const InstanceT& inst) {
        Cnts usage(inst.nb_items(), 0);
        inst.enum_up([&] (VertexID r) {
            inst.enum_hyperarcs(r, [&] (const Hyperarc& e) {
                for(int item_id = 0;item_id < inst.nb_items();item_id++) {
                    if(inst.items[item_id].vrt == e.u) usage[item_id]++;
                    if(inst.items[item_id].vrt == e.v) usage[item_id]++;
                }
            });
        });
        return usage;
    }

    template<typename InstanceT>
    Obj solution_value(const InstanceT& inst) {
        Cnts usage = usage_of_items(inst);
        Obj obj = 0;
        for(int item_id = 0;item_id < inst.nb_items();item_id++) {
            obj += inst.items[item_id].obj * std::min(inst.items[item_id].ub, usage[item_id]);
        }
        return obj;
    }

    template<typename InstanceT>
    std::vector<VertexID> compute_order(const InstanceT& inst) {
        std::vector<VertexID> order;

        enum Status {
            NONE, PENDING, DONE
        };

        std::vector<Status> status(inst.nb_vertices(), NONE);

        auto explore_rec = [&] (auto self, VertexID r) -> void {
            if(status[r] == PENDING) {
                std::cout << "cycle detected" << std::endl;
                exit(0);
            }
            if(status[r] == NONE) {
                status[r] = PENDING;
                inst.enum_hyperarcs(r, [&] (Hyperarc hyperarc) {
                    self(self, hyperarc.u);
                    self(self, hyperarc.v);
                });
                order.push_back(r);
                status[r] = DONE;
            }
        };

        explore_rec(explore_rec, inst.root);
        return order;
    }
}

#endif