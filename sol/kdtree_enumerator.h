// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

#ifndef KDTREE_ENUMERATOR_H
#define KDTREE_ENUMERATOR_H

#include "indexed_label_set.h"
#include "potential_filter.h"

namespace HPG {
    struct ExtendedLabel {
        std::vector<uint8_t> cnts;
        Obj obj;
    };

    const size_t KDTREE_LEAF_SIZE = 10;

    struct KDTree {
        KDTree *left = nullptr, *right = nullptr;
        std::vector<uint8_t> lbs, ubs;
        std::vector<ExtendedLabel*> labels;
        Obj best_leaf = 0;
        
        ~KDTree();
    };

    struct KDTreeEnumerator {
        const PotentialFilter& main_filter;

        KDTreeEnumerator(const PotentialFilter& _main_filter);
        KDTree* build(VertexID r, std::vector<ExtendedLabel*>& labels) const;
        void enum_pairs(Hyperarc e, std::vector<uint64_t>& labels_u, std::vector<uint64_t>& labels_v, uint64_t* label, size_t nb_chunks, float lb, auto func) const;
        void compile(VertexID r, std::vector<uint64_t>& labels, size_t nb_chunks) const;
    };
}

namespace HPG {
    KDTree::~KDTree() {
        if(left != nullptr) delete left;
        if(right != nullptr) delete right;
    }

    KDTree* KDTreeEnumerator::build(VertexID r, std::vector<ExtendedLabel*>& labels) const {
        const int nb_items = main_filter.items.size();

        std::vector<uint8_t> lbs(nb_items, 255);
        std::vector<uint8_t> ubs(nb_items, 0);

        for(auto label : labels) {
            for(int item_id = 0;item_id < nb_items;item_id++) {
                lbs[item_id] = std::min(lbs[item_id], label->cnts[item_id]);
                ubs[item_id] = std::max(ubs[item_id], label->cnts[item_id]);
            }
        }

        if(labels.size() < KDTREE_LEAF_SIZE) {
            std::stable_sort(labels.begin(), labels.end(), [] (ExtendedLabel* a, ExtendedLabel* b) {
                return a->obj > b->obj;
            });
            return new KDTree{nullptr, nullptr, lbs, ubs, labels, labels[0]->obj};
        }

        size_t best_split = 0;
        Obj best_range = 0;
        for(size_t iItem = 0;iItem < main_filter.items.size();iItem++) {
            Obj item_range = (ubs[iItem] - lbs[iItem]) * main_filter.filter_coeff(r, iItem);
            if(item_range > best_range) {
                best_range = item_range;
                best_split = iItem;
            }
        }

        std::stable_sort(labels.begin(), labels.end(), [&best_split] (ExtendedLabel* a, ExtendedLabel* b) {
            return a->cnts[best_split] < b->cnts[best_split];
        });

        size_t median = labels.size() / 2;
        std::vector<ExtendedLabel*> left_labels(labels.begin(), labels.begin() + median);
        std::vector<ExtendedLabel*> right_labels(labels.begin() + median, labels.end());

        KDTree* tree = new KDTree{build(r, left_labels), build(r, right_labels), lbs, ubs};
        tree->best_leaf = std::max(tree->left->best_leaf, tree->right->best_leaf);
        return tree;
    }

    KDTreeEnumerator::KDTreeEnumerator(const PotentialFilter& _main_filter) : main_filter(_main_filter) {
    }

    void KDTreeEnumerator::enum_pairs(Hyperarc e, std::vector<uint64_t>& labels_u, std::vector<uint64_t>& labels_v, uint64_t* label, size_t nb_chunks, float lb, auto func) const {
        if(labels_u.empty() || labels_v.empty()) return;
        
        const int nb_items = main_filter.items.size();

        Cnts ubs(nb_items);
        float coeffs[nb_items];

        for(int item_id = 0;item_id < nb_items;item_id++) {
            ubs[item_id] = main_filter.items[item_id].ub;
            coeffs[item_id] = main_filter.filter_coeff(e.r, item_id);
        }

        auto explore_rec = [&](auto self, ExtendedLabel* label, KDTree* node) -> void {
            Obj max_obj = 0;
            Obj other_bound = node->best_leaf + label->obj;

            for(int item_id = 0;item_id < nb_items;item_id++) {
                max_obj += std::min<Cnt>(label->cnts[item_id] + node->ubs[item_id], main_filter.items[item_id].ub) * coeffs[item_id];
                other_bound -= std::max<int64_t>(0, label->cnts[item_id] + node->lbs[item_id] - (int64_t)main_filter.items[item_id].ub) * coeffs[item_id];
            }

            if(std::min(other_bound, max_obj) + main_filter.constant(e.r) < lb) return;

            if(node->left == nullptr) {
                for(ExtendedLabel* other : node->labels) {
                    if(label->obj + other->obj + main_filter.constant(e.r) < lb) return;

                    Obj score = main_filter.constant(e.r);
                    for(int item_id = 0;item_id < nb_items;item_id++) {
                        score += std::min<Cnt>(label->cnts[item_id] + other->cnts[item_id], main_filter.items[item_id].ub) * coeffs[item_id];
                    }

                    if(score >= lb) {
                        func(other->cnts, label->cnts);
                    }
                }
            } else {
                self(self, label, node->left);
                self(self, label, node->right);
            }
        };
        
        ExtendedLabel best_u;
        best_u.cnts = expand_label(ubs, &labels_u[0], nb_chunks);
        best_u.obj = main_filter.linear(e.u, best_u.cnts);

        ExtendedLabel best_v;
        best_v.cnts = expand_label(ubs, &labels_v[0], nb_chunks);
        best_v.obj = main_filter.linear(e.v, best_v.cnts);

        float constant = main_filter.rdp[e.r] + main_filter.hyperarc_obj(e);

        std::vector<ExtendedLabel> elabels_u;
        for(size_t uLabel = 0;uLabel < labels_u.size();uLabel += nb_chunks) {
            ExtendedLabel uElabel;
            uElabel.cnts = expand_label(ubs, &labels_u[uLabel], nb_chunks);
            uElabel.obj = main_filter.linear(e.r, uElabel.cnts);

            if(main_filter.linear(e.u, uElabel.cnts) + best_v.obj + constant < lb) {
                break;
            }

            elabels_u.push_back(uElabel);
        }

        if(elabels_u.empty()) return;

        std::vector<ExtendedLabel*> ptrs;
        for(size_t iLabel = 0;iLabel < elabels_u.size();iLabel++) {
            ptrs.push_back(&elabels_u[iLabel]);
        }

        auto tree = build(e.r, ptrs);

        for(size_t vLabel = 0;vLabel < labels_v.size();vLabel += nb_chunks) {
            ExtendedLabel vElabel;
            vElabel.cnts = expand_label(ubs, &labels_v[vLabel], nb_chunks);
            vElabel.obj = main_filter.linear(e.r, vElabel.cnts);

            if(main_filter.linear(e.v, vElabel.cnts) + best_u.obj + constant < lb) {
                break;
            }

            explore_rec(explore_rec, &vElabel, tree);
        }

        delete tree;
    }

    void KDTreeEnumerator::compile(VertexID r, std::vector<uint64_t>& labels, size_t nb_chunks) const {
        Cnts ubs(main_filter.items.size());
        for(int item_id = 0;item_id < (int)main_filter.items.size();item_id++) {
            ubs[item_id] = main_filter.items[item_id].ub;
        }

        std::vector<std::pair<Obj, uint64_t*>> obj_labels;
        for(int label_id = 0;label_id < (int)labels.size();label_id += nb_chunks) {
            obj_labels.push_back({-main_filter.linear(r, expand_label(ubs, &labels[label_id], nb_chunks)), &labels[label_id]});
        }

        std::stable_sort(obj_labels.begin(), obj_labels.end());

        std::vector<uint64_t> new_labels;
        for(int label_id = 0;label_id < obj_labels.size();label_id++) {
            new_labels.insert(new_labels.end(), obj_labels[label_id].second, obj_labels[label_id].second + nb_chunks);
        }
        labels = new_labels;
    }
}

#endif