// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains the data structure for a set of labels, encoded as a bitset.
*/

#ifndef INDEXED_LABEL_SET_H
#define INDEXED_LABEL_SET_H

#include "label.h"

namespace HPG {
    const int MAX_BUFFER_SIZE = 256;

    struct IndexedLabelSet {
    private:
        // number of 64-bit integers needed to encode a label 
        int nb_chunks;

        // for each bit, an array containing all the labels with the bit set, stored contiguously
        std::vector<std::vector<uint64_t>> index;

        std::vector<std::vector<uint64_t>> buffers;

        void real_add(uint64_t label[]);
        void clear_buffer(int bit);
    
    public:
        // array containing all the labels stored contiguously
        std::vector<uint64_t> labels;

        // create the data structure
        void init(int _nb_chunks);
        
        // add a label to the data structure
        void add_label(uint64_t label[]);

        // free the data structures, only keeps the labels
        // must be called when all the labels have been added
        void close();
    };
}

namespace HPG {
    void IndexedLabelSet::real_add(uint64_t label[]) {
        labels.insert(labels.end(), label, label + nb_chunks);

        for(size_t iChunk = 0;iChunk < nb_chunks;iChunk++) {
            uint64_t bits = label[iChunk];
            while(bits != 0) {
                int iBit = std::countr_zero(bits);
                index[64 * iChunk + iBit].insert(index[64 * iChunk + iBit].end(), label, label + nb_chunks);
                bits &= bits - 1;
            }
        }
    }

    void IndexedLabelSet::clear_buffer(int bit) {
        for(int iLabel = 0;iLabel < (int)index[bit].size();iLabel += nb_chunks) {
            for(int iOther = 0;iOther < (int)buffers[bit].size();) {
                if(compare_labels(&buffers[bit][iOther], &index[bit][iLabel], nb_chunks)) {
                    std::copy(buffers[bit].end() - nb_chunks, buffers[bit].end(), buffers[bit].begin() + iOther);
                    buffers[bit].resize(buffers[bit].size() - nb_chunks);
                } else {
                    iOther += nb_chunks;
                }
            }
        }

        for(size_t iOther = 0;iOther < buffers[bit].size();iOther += nb_chunks) {
            real_add(&buffers[bit][iOther]);
        }
        buffers[bit].clear();
    }

    void IndexedLabelSet::init(int _nb_chunks) {
        nb_chunks = _nb_chunks;
        buffers.resize(64 * nb_chunks);
        index.resize(64 * nb_chunks);
    }

    void IndexedLabelSet::add_label(uint64_t label[]) {
        size_t best_bit = 0, best_size = INT_MAX;

        for(size_t iChunk = 0;iChunk < nb_chunks;iChunk++) {
            uint64_t bits = label[iChunk];
            while(bits != 0) {
                int iBit = std::countr_zero(bits);

                size_t sz = index[64 * iChunk + iBit].size() + buffers[64 * iChunk + iBit].size();
                if(sz < best_size) {
                    best_size = sz;
                    best_bit = 64 * iChunk + iBit;
                }

                bits &= (bits - 1);
            }
        }

        if(best_size != INT_MAX) {
            for(size_t iLabel = 0;iLabel < buffers[best_bit].size();iLabel += nb_chunks) {
                if(compare_labels(label, &buffers[best_bit][iLabel], nb_chunks)) {
                    return;
                }
            }

            buffers[best_bit].insert(buffers[best_bit].end(), label, label + nb_chunks);

            if(buffers[best_bit].size() > nb_chunks * MAX_BUFFER_SIZE) {
                clear_buffer(best_bit);
            }
        } else if(labels.empty()) {
            real_add(label);
        }
    }

    void IndexedLabelSet::close() {
        for(int iBit = 0;iBit < 64 * nb_chunks;iBit++) {
            clear_buffer(iBit);
        }

        buffers.clear(); buffers.shrink_to_fit();
        index.clear(); index.shrink_to_fit();
    }
}

#endif