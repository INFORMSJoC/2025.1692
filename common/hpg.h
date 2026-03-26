// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains the explicit hypergraph data structure.
*/

#ifndef HPG_H
#define HPG_H

#include "types.h"
#include "hpg_algorithms.h"

namespace HPG {
    // the type of an explicit hypergraph
    struct HPG {
        // the array of items that can be produced
        Items items;

        // for each vertex id, hyperarcs having the vertex as target
        std::vector<Hyperarcs> hyperarcs;
    
        // for each vertex, its name
        std::vector<std::string> vertex_names;

        // id of the vertex of the initial object
        VertexID root = -1;

        // topological ordering of the vertices
        std::vector<VertexID> order;

        // create an empty hypergraph
        HPG();

        // create an hypergraph from a hpg file
        HPG(const std::string& path);
        
        // create a new vertex and return the id
        VertexID new_vertex(const std::string& name);

        // create a new item
        void new_item(const Item& item);

        // create a new hyperarc
        void new_hyperarc(const Hyperarc& hyperarc);

        // return the number of vertices
        size_t nb_vertices() const;

        // return the number of items
        size_t nb_items() const;

        // call the function func for each hyperarc having r for target
        void enum_hyperarcs(VertexID r, auto func) const;

        // call the function func for each vertex id in the bottom-top topological ordering
        void enum_up(auto func) const;

        // call the function func for each vertex id in the top-bottom topological ordering
        void enum_down(auto func) const;
    };

    // return an hypergraph with the same root, vertices and items, but with no hyperarcs
    template<typename InstanceT>
    HPG remove_hyperarcs(const InstanceT& inst);
}

namespace HPG {
    HPG::HPG() {}

    HPG::HPG(const std::string& path) {
        std::ifstream fin;
        fin.open(path);
        if(!fin.is_open()) {
            std::cout << "unable to open " << path << std::endl;
            exit(0);
        }

        auto error = [&] () {
            std::cout << "error while parsing " << path << std::endl;
            exit(0);
        };

        enum Section {
            COMMENT, SYMBOLS, EDGES, ROOT, ITEMS
        };

        Section section = COMMENT;
        VertexID current_r;

        while(true) {
            std::string line;
            getline(fin, line);

            if(fin.fail()) {
                if(section == ITEMS) break;
                else error();
            }

            if(line == "<symbols>") {
                if(section == COMMENT) section = SYMBOLS;
                else error();
            } else if(line == "<edges>") {
                if(section == SYMBOLS) section = EDGES;
                else error();
            } else if(line == "<root>") {
                if(section == EDGES) section = ROOT;
                else error();
            } else if(line == "<items>") {
                if(section == ROOT) section = ITEMS;
                else error();
            } else {
                std::stringstream ss;
                ss << line;

                if(section == COMMENT) {
                } else if(section == ROOT) {
                    ss >> root;
                    if(ss.fail()) error();
                } else if(section == EDGES) {
                    VertexID first_id;
                    ss >> first_id;

                    if(ss.fail()) error();

                    VertexID second_id;
                    ss >> second_id;

                    if(ss.fail()) current_r = first_id;
                    else new_hyperarc({current_r, first_id, second_id});
                } else if(section == ITEMS) {
                    Item item;
                    ss >> item.vrt >> item.obj >> item.ub;
                    if(ss.fail()) error();
                    new_item(item);
                } else if(section == SYMBOLS) {
                    std::string name;
                    ss >> name;
                    if(ss.fail()) error();
                    new_vertex(name);
                }
            }
        }

        if(root == -1) error();

        fin.close();
        order = compute_order(*this);
    }

    VertexID HPG::new_vertex(const std::string& name) {
        vertex_names.push_back(name);
        hyperarcs.push_back({});
        return vertex_names.size() - 1;
    }

    inline void HPG::new_item(const Item& item) {
        items.push_back(item);
    }

    void HPG::new_hyperarc(const Hyperarc& hyperarc) {
        hyperarcs[hyperarc.r].push_back(hyperarc);
    }

    inline size_t HPG::nb_vertices() const {
        return vertex_names.size();
    }

    inline size_t HPG::nb_items() const {
        return items.size();
    }

    void HPG::enum_hyperarcs(VertexID r, auto func) const {
        for(const auto& hyperarc : hyperarcs[r]) {
            func(hyperarc);
        }
    }

    void HPG::enum_up(auto func) const {
        for(const auto& vertex_id : order) {
            func(vertex_id);
        }
    }

    void HPG::enum_down(auto func) const {
        for (auto it = order.rbegin();it != order.rend();it++) {
            func(*it);
        }
    }

    template<typename InstanceT>
    HPG remove_hyperarcs(const InstanceT& inst) {
        HPG hpg;
        hpg.root = inst.root;
        hpg.items = inst.items;
        hpg.vertex_names = inst.vertex_names;
        hpg.hyperarcs.resize(hpg.vertex_names.size());
        hpg.order = inst.order;
        return hpg;
    }
}

#endif