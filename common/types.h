// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains all the needed STL includes and basic types.
*/

#ifndef TYPES_H
#define TYPES_H

#include "stdint.h"
#include <vector>
#include <bitset>
#include <algorithm>
#include <climits>
#include <cmath>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <functional>
#include <deque>
#include <numeric>
#include <thread>
#include <chrono>
#include <cassert>

namespace HPG {
    // type of a vertex id
    using VertexID = int32_t;

    // type of a profit
    using Obj = float;

    // type of an array of profits
    using Objs = std::vector<float>;

    // type of an item upper bound
    using Cnt = int32_t;

    // type of an array of item upper bounds
    using Cnts = std::vector<Cnt>;

    // type of an hyperarc
    struct Hyperarc {
        // the target of the hyperarc
        VertexID r;

        // the two sources of the hyperarc
        VertexID u, v;
    };

    // type of an array of hyperarcs
    using Hyperarcs = std::vector<Hyperarc>;

    // comparison operator for hyperarcs
    bool operator < (const Hyperarc& a, const Hyperarc& b);

    // type of an item
    struct Item {
        // the vertex id where the item can be produced
        VertexID vrt;

        // the profit a single item
        Obj obj;

        // the number of times the item can be produced    
        Cnt ub;
    };

    // type of an array of items
    using Items = std::vector<Item>;
}

namespace HPG {
    // we compare target and source ids lexicographically
    bool operator < (const Hyperarc& a, const Hyperarc& b) {
        if(a.r != b.r) return a.r < b.r;
        if(a.u != b.u) return a.u < b.u;
        return a.v < b.v;
    }
}

#endif
