// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

#ifndef PERMUTATION_H
#define PERMUTATION_H

#include "../common/hpg_algorithms.h"

namespace HPG {
    using Permutation = std::vector<int>;

    template<typename T>
    void permute(const Permutation& permutation, const std::vector<T>& input, std::vector<T>& output);

    std::vector<Permutation> all_permutations(size_t dim);
}

namespace HPG {
    template<typename T>
    inline void permute(const Permutation& permutation, const std::vector<T>& input, std::vector<T>& output) {
        for(int id = 0;id < (int)input.size();id++) {
            output[permutation[id]] = input[id];
        }
    }

    std::vector<Permutation> all_permutations(size_t dim) {
        std::vector<Permutation> permutations;
        Permutation permutation(dim, 0);
        std::iota(permutation.begin(), permutation.end(), 0);
        do {
            permutations.push_back(permutation);
        } while(std::next_permutation(permutation.begin(), permutation.end()));
        return permutations;
    }
}

#endif