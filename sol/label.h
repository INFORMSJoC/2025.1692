// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains the algorithms to manipulate labels.
    Labels can be encoded in two different ways : 
    - as an array, in a uint8_t arrat  : for each item the number of times it is in the label 
    - as a bitset, in a uint64_t array : for each possible item instance, if it is present or not in the label 
*/

#ifndef LABEL_H
#define LABEL_H

#include "../common/hpg.h"

namespace HPG {
    size_t nb_chunks_of(const Items& items);
    void empty_label(uint64_t label[], size_t nb_chunks);
    void item_label(size_t iItem, const Items& items, uint64_t label[], size_t nb_chunks);
    bool compare_labels(uint64_t a[], uint64_t b[], size_t nb_chunks);
    std::vector<uint8_t> expand_label(const Items& items, uint64_t label[], size_t nb_chunks);
    void merge_labels(const Items& items, const std::vector<uint8_t>& a, const std::vector<uint8_t>& b, uint64_t label[], size_t nb_chunks);
    Obj label_obj(const Items& items, const std::vector<uint8_t>& label);
}

namespace HPG {
    size_t nb_chunks_of(const Items& items) {
        size_t nb_items = 0;
        for(const Item& item : items) {
            nb_items += item.ub;
        }
        return (nb_items + 63) / 64;
    }

    void empty_label(uint64_t label[], size_t nb_chunks) {
        std::fill(label, label + nb_chunks, 0);
    }

    void item_label(size_t iItem, const Items& items, uint64_t label[], size_t nb_chunks) {
        std::fill(label, label + nb_chunks, 0);
        
        size_t shift = 0;
        for(size_t iBefore = 0;iBefore < iItem;iBefore++) {
            shift += items[iBefore].ub;
        }

        label[shift / 64] |= (uint64_t)1 << (shift % 64);
    }

    inline bool compare_labels(uint64_t a[], uint64_t b[], size_t nb_chunks) {
        for(size_t iChunk = 0;iChunk < nb_chunks;iChunk++) {
            if((a[iChunk] | b[iChunk]) != b[iChunk]) {
                return false;
            }
        }
        return true;
    }

    std::vector<uint8_t> expand_label(const Cnts& ubs, uint64_t label[], size_t nb_chunks) {
        std::vector<uint8_t> cnts(ubs.size(), 0);

        size_t item = 0, base = 0;
        for(size_t iChunk = 0;iChunk < nb_chunks;iChunk++) {
            uint64_t bits = label[iChunk];

            while(bits != 0) {
                int iBit = std::countr_zero(bits);

                size_t bit_id = 64 * iChunk + iBit;
                while(bit_id - base >= ubs[item]) {
                    base += ubs[item];
                    item++;
                }
                cnts[item] += 1;

                bits &= (bits - 1);
            }
        }

        return cnts;
    }

    void merge_labels(const Items& items, const std::vector<uint8_t>& a, const std::vector<uint8_t>& b, uint64_t label[], size_t nb_chunks) {
        std::fill(label, label + nb_chunks, 0);

        size_t base = 0;
        for(size_t iItem = 0;iItem < items.size();iItem++) {
            auto num = std::min<uint8_t>(a[iItem] + b[iItem], items[iItem].ub);
            for(size_t iBit = base;iBit < base + num;iBit++) {
                label[iBit / 64] |= (uint64_t)1 << (iBit % 64);
            }
            base += items[iItem].ub;
        }
    }

    Obj label_obj(const Items& items, const std::vector<uint8_t>& label) {
        Obj obj = 0;
        for(size_t iItem = 0;iItem < items.size();iItem++) {
            obj += items[iItem].obj * label[iItem];
        }
        return obj;
    }
}

#endif