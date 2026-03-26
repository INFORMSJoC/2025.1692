// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains all the functions that are useful for parsing
    a geometric instance file, and transforming it into an explicit decision hypergraph.
*/

#ifndef GEOMETRIC_INSTANCE_FILE_H
#define GEOMETRIC_INSTANCE_FILE_H

#include "encoder_arguments.h"

namespace HPG {
    using Box = std::vector<int>;

    struct GeometricInstanceFile {
        // dimensions of the initial plate
        Box initial_box;

        // dimensions of the items
        std::vector<Box> item_boxes;

        // profits of the items
        std::vector<Obj> item_objs;

        // maximum number of times each item can be sold
        std::vector<Cnt> item_ubs;

        // encoder arguments
        EncoderArguments arguments;

        // implicit hypergraph arguments
        std::vector<std::vector<bool>> is_restricted_for;
        std::vector<bool> is_useful_box;
        std::vector<int> loose_axes;
        std::vector<std::vector<int>> normalizations;

        std::vector<std::vector<VertexID>> vertices;

        std::vector<std::string> vertex_names;
        std::vector<Box> boxes;
        std::vector<int> states;

        Items items;
        VertexID root = -1;

        std::vector<VertexID> order;

        // parse a geometric instance file
        GeometricInstanceFile(const EncoderArguments& _arguments);

        // transform the geometric instance into a decision hypergraph
        void hpg();

        // convert a box into a string representation
        std::string str(const Box& box);
        
        // convert a vertex into a string representation
        std::string str(const Box& box, int state);

        // call func for each possible box
        void iter_boxes(auto func);

        // convert a box_id to the box
        void box_of(Box& box, size_t box_id);

        // return the box_id corresponding to a box
        int id_of(const Box& box) const;

        // check if a box is contained in an other box
        bool is_contained(const Box& box1, const Box& box2) const;

        // return the measure of a box
        int measure_of(const Box& box) const;

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
}

namespace HPG {
    GeometricInstanceFile::GeometricInstanceFile(const EncoderArguments& _arguments) : arguments(_arguments) {
        // open the input file
        std::ifstream fin;
        fin.open(arguments.input);
        if(!fin.is_open()) {
            std::cerr << "unable to open " << arguments.input << std::endl;
            exit(0);
        }

        auto error = [&] () {
            std::cerr << "error while parsing " << arguments.input << std::endl;
            exit(0);
        };

        // read the dimensions of the initial plate
        initial_box.resize(arguments.dim);

        int initial_measure = 1;
        for(int axis = 0;axis < arguments.dim;axis++) {
            int coord;
            fin >> coord;
            if(fin.fail()) error();

            initial_box[axis] = coord / arguments.reduction_ratio[axis];
            initial_measure *= initial_box[axis];
            assert(initial_box[axis] > 0);
        }

        // read the number of items
        size_t nb_items;
        fin >> nb_items;
        if(fin.fail()) error();

        // read the data of the items : dimensions, profit, sales bound
        item_boxes.resize(nb_items); item_objs.resize(nb_items); item_ubs.resize(nb_items);
        for(size_t item_id = 0;item_id < nb_items;item_id++) {
            item_boxes[item_id].resize(arguments.dim);

            int measure = 1, reduced_measure = 1;
            for(int axis = 0;axis < arguments.dim;axis++) {
                int coord;
                fin >> coord;
                if(fin.fail()) error();

                if(arguments.roundup) item_boxes[item_id][axis] = (coord + arguments.reduction_ratio[axis] - 1) / arguments.reduction_ratio[axis];
                else item_boxes[item_id][axis] = coord / arguments.reduction_ratio[axis];

                assert(item_boxes[item_id][axis] > 0);
                measure *= coord;
                reduced_measure *= item_boxes[item_id][axis];
            }

            fin >> item_objs[item_id] >> item_ubs[item_id];

            item_ubs[item_id] = std::min(item_ubs[item_id], initial_measure / reduced_measure);

            // if the --measure option is activated, update the profits with the measure (length, area, volume, ...)
            if(arguments.measure) item_objs[item_id] = measure;
            if(fin.fail()) error();
        }
    
        fin.close();
    }

    std::string GeometricInstanceFile::str(const Box& box) {
        std::string name = std::to_string(box[0]);
        for(size_t dim = 1;dim < box.size();dim++) {
            name += "," + std::to_string(box[dim]);
        }
        return name;
    }

    std::string GeometricInstanceFile::str(const Box& box, int state) {
        return str(box) + "/" + std::to_string(state);
    }

    void GeometricInstanceFile::iter_boxes(auto func) {
        Box box(arguments.dim, 0);

        auto iter_rec = [&](auto self, int axis) -> void {
            if(axis == arguments.dim) {
                func(box);
                return;
            }
            
            for(int length = 1;length <= initial_box[axis];length++) {
                box[axis] = length;
                self(self, axis + 1);
            }
        };

        iter_rec(iter_rec, 0);
    };

    inline void GeometricInstanceFile::box_of(Box& box, size_t box_id) {
        for(int axis = 0;axis < initial_box.size();axis++) {
            box[axis] = box_id % (1 + initial_box[axis]);
            box_id /= (1 + initial_box[axis]);
        }
    }

    inline int GeometricInstanceFile::id_of(const Box& box) const {
        int base = 1, state = 0;
        for(int axis = 0;axis < arguments.dim;axis++) {
            state += base * box[axis];
            base *= (1 + initial_box[axis]);
        }
        return state;
    }

    inline bool GeometricInstanceFile::is_contained(const Box& box1, const Box& box2) const {
        for(int axis = 0;axis < arguments.dim;axis++) {
            if(box1[axis] > box2[axis]) return false;
        }
        return true;
    }

    inline int GeometricInstanceFile::measure_of(const Box& box) const {
        int measure = 1;
        for(int length : box) measure *= length;
        return measure;
    }

    void GeometricInstanceFile::hpg() {
        // compute the total number of boxes (including zeroes)
        int nb_boxes = 1;
        for(int length : initial_box) nb_boxes *= (1 + length);

        // for each possible box, compute if it can be removed, because it is not
        // the minimal size of some cutting pattern which satisfy the sales bounds.
        // we start by marking all the boxes as useful.
        is_useful_box = std::vector<bool>(nb_boxes, true);

        // we create a global buffer to avoid reallocating the memory each iteration
        Box sum_buffer(arguments.dim);
        Box box(arguments.dim);
        Box normalized_box(arguments.dim);
        Box reduced_box(arguments.dim);
        Box rotated_item(arguments.dim);

        // reduce item ubs based on domination arguments
        if(std::count(arguments.trim.begin(), arguments.trim.end(), false) == 0) {
            int delta = 0;
            for(int item_id = 0;item_id < (int)item_boxes.size();item_id++) {
                int remaining_measure = measure_of(initial_box);
                for(int other_item_id = 0;other_item_id < (int)item_boxes.size();other_item_id++) {
                    if(item_id == other_item_id || item_objs[other_item_id] <= item_objs[item_id]) continue;
                    
                    bool is_dominated = false;
                    for(const auto& permutation : arguments.rotations) {
                        permute(permutation, item_boxes[other_item_id], rotated_item);
                        is_dominated |= is_contained(rotated_item, item_boxes[item_id]);
                    }

                    if(is_dominated) {
                        remaining_measure = std::max(0, remaining_measure - item_ubs[other_item_id] * measure_of(item_boxes[other_item_id]));
                    }
                }

                item_ubs[item_id] = std::min(item_ubs[item_id], remaining_measure / measure_of(item_boxes[item_id]));
            }
        }

        // for each axis, we use dynamic programming to compute the smallest dimension along that axis
        // of an item box that must appear in a sum giving a box.
        for(int measured_axis = 0;measured_axis < arguments.dim;measured_axis++) {
            std::vector<int> smallest_for(nb_boxes, INT_MAX);
            smallest_for[0] = 0;

            // for each copy of each item, we update the DP table (smallest_for)
            for(size_t item_id = 0;item_id < item_boxes.size();item_id++) {
                for(size_t _ = 0;_ < item_ubs[item_id];_++) {
                    for(int box_id = nb_boxes - 1;box_id >= 0;box_id--) {
                        box_of(box, box_id);

                        // each copy can be rotated
                        for(const auto& permutation : arguments.rotations) {
                            permute(permutation, item_boxes[item_id], rotated_item);

                            // each copy can contribute to every non-empty set of axis in the sum
                            for(int mask = 1;mask < (1 << arguments.dim);mask++) {
                                // compute the sum
                                for(int axis = 0;axis < arguments.dim;axis++) {
                                    sum_buffer[axis] = box[axis] + ((mask & (1 << axis)) != 0) * rotated_item[axis];
                                }

                                // update the DP table if the obtained sum is contained in the initial plate
                                if(is_contained(sum_buffer, initial_box)) {
                                    int sum_id = id_of(sum_buffer);
                                    smallest_for[sum_id] = std::min(smallest_for[sum_id], std::max(smallest_for[box_id], rotated_item[measured_axis]));
                                }
                            }     
                        }
                    }
                }
            }
            
            // mark as useless boxes for which no sum can be written
            for(int box_id = 0;box_id < nb_boxes;box_id++) {
                box_of(box, box_id);
                is_useful_box[box_id] = is_useful_box[box_id] && (smallest_for[box_id] <= box[measured_axis]);
            }
        }

        // compute for each box, the loose axes : the axes parallel to which a solution can always be reduced
        loose_axes = std::vector<int>(nb_boxes, 0);
        iter_boxes([&] (Box& box) {
            int box_id = id_of(box);
            // if a box can be expressed as a sum, it has no loose axis
            if(!is_useful_box[box_id]) {
                // otherwise, an axis is loose if it is loose for every reduced box
                loose_axes[box_id] = (1 << arguments.dim) - 1;
                for(int axis = 0;axis < arguments.dim;axis++) {
                    if(box[axis] >= 2) {
                        box[axis] -= 1;
                        int reduced_box_id = id_of(box);
                        box[axis] += 1;
                        loose_axes[box_id] &= loose_axes[reduced_box_id] | (1 << axis);
                    }
                }
            }
        });

        // for each axis and each dimension, the normalized dimension
        for(int axis = 0;axis < arguments.dim;axis++) {
            std::vector<int> normalized(1 + initial_box[axis], 0);
            normalized[initial_box[axis]] = initial_box[axis];
            iter_boxes([&] (Box& box) {
                int box_id = id_of(box);
                if(is_useful_box[box_id]) {
                    normalized[box[axis]] = box[axis];
                }
            });
            for(int dim = 1;dim <= initial_box[axis];dim++) {
                normalized[dim] = std::max(normalized[dim], normalized[dim - 1]);
            }
            normalizations.push_back(normalized);
        }

        // compute restricted cut lengths
        for(int axis = 0;axis < arguments.dim;axis++) {
            std::vector<bool> is_restricted(1 + initial_box[axis], false);
            for(size_t item_id = 0;item_id < item_boxes.size();item_id++) {
                for(const auto& permutation : arguments.rotations) {
                    permute(permutation, item_boxes[item_id], rotated_item);

                    if(is_contained(rotated_item, initial_box)) {
                        is_restricted[rotated_item[axis]] = true;
                    }
                }
            }
            is_restricted_for.push_back(is_restricted);
        }

        // create the special none vertex
        vertex_names.push_back("none");
        boxes.push_back({});
        states.push_back(0);

        // create the vertices of each item and its rotations
        for(size_t item_id = 0;item_id < item_boxes.size();item_id++) {
            std::string name = "i" + std::to_string(item_id);
            items.push_back(Item{(int)vertex_names.size(), item_objs[item_id], item_ubs[item_id]});
            vertex_names.push_back(name);
            boxes.push_back({});
            states.push_back(0);
        }

        // create the other vertices
        vertices = std::vector<std::vector<VertexID>>(arguments.stages.states.size(), std::vector<VertexID>(nb_boxes, -1));
        for(size_t state_id = 0;state_id < arguments.stages.states.size();state_id++) {
            iter_boxes([&] (Box& box) {
                int box_id = id_of(box);
                bool is_normalized = true;
                for(int axis = 0;axis < arguments.dim;axis++) {
                    is_normalized &= normalizations[axis][box[axis]] == box[axis];
                }

                if(is_normalized) {
                    vertices[state_id][box_id] = vertex_names.size();
                    vertex_names.push_back(str(box, state_id));
                    boxes.push_back(box);
                    states.push_back(state_id);
                }
            });
        }

        // create the root
        root = vertices[arguments.stages.start][id_of(initial_box)];

        // update the topological order of the vertices
        order = compute_order(*this);
    }

    void GeometricInstanceFile::enum_hyperarcs(VertexID r, auto func) const {
        if(vertex_names[r] == "none" || vertex_names[r][0] == 'i') return;

        // hyperarcs boxes -> items
        Box rotated_item(arguments.dim);
        for(size_t item_id = 0;item_id < item_boxes.size();item_id++) {
            for(const auto& permutation : arguments.rotations) {
                permute(permutation, item_boxes[item_id], rotated_item);

                if(is_contained(rotated_item, initial_box)) {
                    VertexID dst = vertices[arguments.stages.end][id_of(rotated_item)];

                    if(dst == r) {
                        func(Hyperarc{dst, items[item_id].vrt, 0});
                    }
                }
            }
        }

        int state_id = states[r];
        const StageGrammarState& state = arguments.stages.states[state_id];

        Box box = boxes[r];
        int box_id = id_of(box);
        Box box_u(arguments.dim), box_v(arguments.dim);

        for(int axis : state.axes) {
            if(is_useful_box[box_id]) {
                // if it is useful, we add every hyperarc
                for(size_t length = 1;length <= box[axis] / 2;length++) {
                    box_u = box; box_u[axis] = length;
                    int box_u_id = id_of(box_u);

                    box_v = box; box_v[axis] = box[axis] - length;
                    int box_v_id = id_of(box_v);

                    // if one of the obtained boxes has the cutting axis loose, we can ignore the hyperarc
                    if((loose_axes[box_u_id] | loose_axes[box_v_id]) & (1 << axis)) continue;

                    // if one axis is loose for both obtained boxes, we can ignore the hyperarc
                    if((loose_axes[box_u_id] & loose_axes[box_v_id]) != 0) continue;
                    
                    // TODO : fix first stage optimisation
                    // if(state.first_stage_optim && box[axis] != initial_box[axis] && box_u[axis] > initial_box[axis] - box[axis]) continue;

                    // last stage optimisation
                    if(state.last_stage_optim && !is_restricted_for[axis][box_u[axis]]) continue;

                    func(Hyperarc{
                        r, 
                        vertices[state_id][box_u_id], 
                        vertices[state_id][box_v_id]
                    });
                }
            } 

            int tmp = box[axis];
            box[axis] = normalizations[axis][tmp - 1];
            
            if(box[axis] != 0) {
                func({
                    r,
                    vertices[state_id][id_of(box)],
                    0
                });
            }
            box[axis] = tmp;
        }

        // create the hyperarcs for the transitions between stages
        for(const auto& transition : arguments.stages.transitions) {
            if(transition.input == state_id) {
                func(Hyperarc{
                    r,
                    vertices[transition.output][box_id],
                    0
                });
            }
        }

        // create the hyperarcs for the trimming
        if(state_id == arguments.stages.end) {
            for(int axis = 0;axis < arguments.dim;axis++) {
                int tmp = box[axis];
                box[axis] = normalizations[axis][tmp - 1];

                if(box[axis] != 0 && arguments.trim[axis]) {
                    func(Hyperarc{
                        r, 
                        vertices[arguments.stages.end][id_of(box)], 
                        0
                    });
                }

                box[axis] = tmp;
            }
        }
    }

    inline size_t GeometricInstanceFile::nb_vertices() const {
        return vertex_names.size();
    }

    inline size_t GeometricInstanceFile::nb_items() const {
        return items.size();
    }

    void GeometricInstanceFile::enum_up(auto func) const {
        for(const auto& vertex_id : order) {
            func(vertex_id);
        }
    }

    void GeometricInstanceFile::enum_down(auto func) const {
        for (auto it = order.rbegin();it != order.rend();it++) {
            func(*it);
        }
    }
}

#endif