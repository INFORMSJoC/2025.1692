// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

#ifndef EXACT_H
#define EXACT_H

#include "kdtree_enumerator.h"
#include "lp.h"

#include <iomanip>

namespace HPG {
    template<typename InstanceT, typename FilterT, typename EnumeratorT>
    HPG solve(const InstanceT& inst, const EnumeratorT& enumerator, Obj lb, int max_size, bool& stopping_condition);

    template<typename InstanceT>
    HPG solve(InstanceT& inst);

    template<typename InstanceT>
    HPG solve_no_cuts(InstanceT& inst);
}

namespace HPG {
    template<typename InstanceT, typename EnumeratorT>
    HPG solve(const InstanceT& inst, const EnumeratorT& enumerator, Obj lb, int max_size, bool& stopping_condition) {
        size_t nb_chunks = nb_chunks_of(inst.items);

        uint64_t label_buffer[nb_chunks];
        std::vector<IndexedLabelSet> labels(inst.nb_vertices());

        logger.max_size = 0;
        logger.total_size = 0;

        Cnts ubs(inst.nb_items());
        for(size_t iItem = 0;iItem < inst.nb_items();iItem++) {
            ubs[iItem] = inst.items[iItem].ub;
        }

        inst.enum_up([&] (VertexID r) {
            if(stopping_condition || enumerator.main_filter.vertex_ub(r) < lb) return;

            labels[r].init(nb_chunks);

            inst.enum_hyperarcs(r, [&](Hyperarc e) {
                if(stopping_condition || enumerator.main_filter.hyperarc_ub(e) < lb) return;

                enumerator.enum_pairs(e, labels[e.u].labels, labels[e.v].labels, label_buffer, nb_chunks, lb, 
                    [&] (const std::vector<uint8_t>& cnts_u, const std::vector<uint8_t>& cnts_v) {
                        merge_labels(inst.items, cnts_u, cnts_v, label_buffer, nb_chunks);
                        labels[r].add_label(label_buffer);
                    }
                );
            });

            for(int item_id = 0;item_id < inst.nb_items();item_id++) {
                if(inst.items[item_id].vrt == r) {
                    item_label(item_id, inst.items, label_buffer, nb_chunks);
                    labels[r].add_label(label_buffer);
                }
            }

            if(r == inst.root || enumerator.main_filter.rdp[r] >= lb) {
                empty_label(label_buffer, nb_chunks);
                labels[r].add_label(label_buffer);
            }

            labels[r].close();

            size_t nb_labels = labels[r].labels.size() / nb_chunks;
            logger.max_size = std::max(logger.max_size, nb_labels);
            logger.total_size += nb_labels;

            if(nb_labels > max_size) {
                std::vector<std::pair<Obj, uint64_t*>> obj_labels;
                for(int label_id = 0;label_id < (int)labels[r].labels.size();label_id += nb_chunks) {
                    auto label = expand_label(ubs, &labels[r].labels[label_id], nb_chunks);
                    Obj score = 0;
                    for(int item_id = 0;item_id < inst.items.size();item_id++) {
                        score += label[item_id] * inst.items[item_id].obj;
                    }
                    obj_labels.push_back(std::pair<Obj, uint64_t*>(-score, &labels[r].labels[label_id]));
                }

                std::stable_sort(obj_labels.begin(), obj_labels.end());

                std::vector<uint64_t> new_labels;
                for(int label_id = 0;label_id < max_size;label_id++) {
                    new_labels.insert(new_labels.end(), obj_labels[label_id].second, obj_labels[label_id].second + nb_chunks);
                }
                labels[r].labels = new_labels;
            }

            enumerator.compile(r, labels[r].labels, nb_chunks);
        });

        if(stopping_condition) return HPG();

        Obj best = 0;
        uint64_t* best_label = &labels[inst.root].labels[0];
        for(size_t rLabel = 0;rLabel < labels[inst.root].labels.size();rLabel += nb_chunks) {
            auto mults = expand_label(ubs, &labels[inst.root].labels[rLabel], nb_chunks);

            Obj obj = 0;
            for(size_t iItem = 0;iItem < inst.nb_items();iItem++) {
                obj += inst.items[iItem].obj * mults[iItem];
            }

            if(obj > best) {
                best = obj;
                best_label = &labels[inst.root].labels[rLabel];
            }
        }

        if(best != 0) {
            logger.lb = std::max(logger.lb, best);
        }

        if(best < lb) return HPG();

        HPG sol = remove_hyperarcs(inst);

        auto create_solution_rec = [&](auto self, VertexID r, std::vector<uint8_t> cnts) -> void {
            if(all_of(cnts.begin(), cnts.end(), [&] (uint8_t mult) { return mult == 0; })) return;

            Hyperarc best_edge = {r, r, r};
            std::vector<uint8_t> best_cnts_u;

            inst.enum_hyperarcs(r, [&](Hyperarc e, Obj obj = 0) {
                enumerator.enum_pairs(e, labels[e.u].labels, labels[e.v].labels, label_buffer, nb_chunks, lb, 
                [&] (const std::vector<uint8_t>& cnts_u, const std::vector<uint8_t>& cnts_v) {
                    for(size_t iItem = 0;iItem < inst.nb_items();iItem++) {
                        if(cnts_u[iItem] + cnts_v[iItem] < cnts[iItem]) {
                            return;
                        }
                    }

                    best_edge = e;
                    best_cnts_u = cnts_u;
                });
            });

            if(best_edge.u == r) return;

            std::vector<uint8_t> cnts_u(cnts.size()), cnts_v(cnts.size());
            for(size_t iItem = 0;iItem < inst.nb_items();iItem++) {
                cnts_u[iItem] = std::min(cnts[iItem], best_cnts_u[iItem]);
                cnts_v[iItem] = cnts[iItem] - cnts_u[iItem];
            }

            self(self, best_edge.u, cnts_u);
            self(self, best_edge.v, cnts_v);

            sol.new_hyperarc(best_edge);
        };

        create_solution_rec(create_solution_rec, inst.root, expand_label(ubs, best_label, nb_chunks));

        sol.order = compute_order(sol);
        return sol;
    }

    template<typename InstanceT>
    HPG solve(InstanceT& inst) {
        logger.start();

        std::mutex mtx;
        bool must_stop = false;

        LpThread lp_thread_obj(inst);

        PotentialFilter filter = lp_thread_obj.hyperarc_generation();

        Obj starting_bound = filter.global_ub() + 1e-3;
        logger.ub = logger.lp_ub = starting_bound;

        logger.lp_event();

        bool integral = true;
        HPG lp_sol;

        PotentialFilter current_filter = filter;

        std::thread lp_thread([&] {
            while(!lp_thread_obj.must_stop) {
                int nb_cuts_before = lp_thread_obj.model.cut_sets.size();
                PotentialFilter improved_filter = lp_thread_obj.cut_generation();
                int nb_cuts_after = lp_thread_obj.model.cut_sets.size();

                if(lp_thread_obj.must_stop) break;

                LpSolution frac_sol = lp_thread_obj.model.get_solution();
                bool is_integer = true;
                for(const auto& it : frac_sol.hyperarc_flow) {
                    is_integer &= fabs(it.second - round(it.second)) <= 1e-5;
                }

                if(is_integer) {
                    logger.event("LP INTEGRAL");

                    std::lock_guard _(mtx);
                    lp_sol = remove_hyperarcs(inst);
                    for(const auto& it : frac_sol.hyperarc_flow) {
                        for(int flow = 0;flow < round(it.second);flow++) {
                            lp_sol.new_hyperarc(it.first);
                        }
                    }
                    must_stop = true;
                    break;
                }

                if(nb_cuts_before == nb_cuts_after) {
                    logger.event("LP STOP");
                    break;
                }

                if(lp_thread_obj.must_stop) break;

                std::lock_guard _(mtx);
                filter.copy(improved_filter);
                
                logger.lp_ub = improved_filter.global_ub();
                if(improved_filter.global_ub() <= 0.9995 * starting_bound) {
                    must_stop = true;
                    logger.lp_event(true);
                } else {
                    logger.lp_event();
                }
            }
        });

        Obj ub = filter.global_ub() - 1e-2;

        int heuristic_size = 3;
        while(true) {
            HPG heuristic_sol = solve(inst, KDTreeEnumerator(current_filter), ub, heuristic_size, must_stop);
            logger.lb = solution_value(heuristic_sol);
            logger.event("HEURISTIC " + std::to_string(heuristic_size));

            if(logger.lb >= logger.ub - 1e-2) {
                lp_thread_obj.must_stop = true;
                lp_thread.join();
                return heuristic_sol;
            }
            
            if(logger.max_size <= heuristic_size) break;

            heuristic_size *= 3;
        }

        while(true) {
            HPG sol;

            do {
                mtx.lock();
                if(lp_sol.nb_vertices() != 0) {
                    sol = lp_sol;
                    break;
                }

                must_stop = false;
                current_filter.copy(filter);
                starting_bound = current_filter.global_ub();

                logger.ub = std::min(logger.ub, starting_bound);
                ub = std::min(ub, starting_bound);
                mtx.unlock();

                sol = solve(inst, KDTreeEnumerator(current_filter), ub - 1e-2, INT_MAX, must_stop);
            } while(must_stop);
            
            if(sol.nb_vertices() != 0) {
                logger.lb = logger.ub = solution_value(sol);
            } else {
                logger.ub = std::min(logger.ub, ub);
            }

            logger.label_event();

            if(sol.nb_vertices() != 0) {
                lp_thread_obj.must_stop = true;
                lp_thread.join();
                return sol;
            }
            
            ub = std::max(logger.lb - 1e-2, 0.999 * ub);
        }
    }

    template<typename InstanceT>
    HPG solve_no_cuts(InstanceT& inst) {
        logger.start();
        bool must_stop = false;

        LpThread<InstanceT>* lp_thread_obj = new LpThread(inst);
        PotentialFilter current_filter = lp_thread_obj->hyperarc_generation();
        
        Obj starting_bound = current_filter.global_ub();
        logger.ub = logger.lp_ub = starting_bound;
        logger.lp_event();
        LpSolution frac_sol = lp_thread_obj->model.get_solution();
        delete lp_thread_obj;

        bool is_integer = true;
        for(const auto& it : frac_sol.hyperarc_flow) {
            is_integer &= fabs(it.second - round(it.second)) <= 1e-5;
        }

        if(is_integer) {
            logger.event("LP INTEGRAL");

            HPG lp_sol = remove_hyperarcs(inst);
            for(const auto& it : frac_sol.hyperarc_flow) {
                for(int flow = 0;flow < round(it.second);flow++) {
                    lp_sol.new_hyperarc(it.first);
                }
            }
            return lp_sol;
        }

        Obj ub = current_filter.global_ub() - 1e-2;

        int heuristic_size = 3;
        while(true) {
            HPG heuristic_sol = solve(inst, KDTreeEnumerator(current_filter), ub, heuristic_size, must_stop);
            logger.lb = solution_value(heuristic_sol);
            logger.event("HEURISTIC " + std::to_string(heuristic_size));

            if(logger.lb >= logger.ub - 1e-2) {
                return heuristic_sol;
            }
            
            if(logger.max_size <= heuristic_size) break;

            heuristic_size *= 3;
        }

        while(true) {
            HPG sol = solve(inst, KDTreeEnumerator(current_filter), ub - 1e-2, INT_MAX, must_stop);
            
            if(sol.nb_vertices() != 0) {
                logger.lb = logger.ub = solution_value(sol);
            } else {
                logger.ub = std::min(logger.ub, ub);
            }

            logger.label_event();

            if(sol.nb_vertices() != 0) {
                return sol;
            }
            
            ub = std::max(logger.lb - 1e-2, 0.999 * ub);
        }
    }
}

#endif