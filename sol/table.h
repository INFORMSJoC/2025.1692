// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains a simple implementation of a dynamically-sized 2D table 
    using a contiguous segment in memory.
*/

#ifndef TABLE_H
#define TABLE_H

#include "../common/hpg_algorithms.h"

namespace HPG {
    template<typename ElemT>
    struct Table {
        // number of rows
        int rows;
        
        // number of columns
        int cols;

        // linearisation of the 2D table
        std::vector<ElemT> table_linearisation;

        // initialise a 2D table with _rows rows and _cols columns, where every cell contains _default
        Table(int _rows = 0, int _cols = 0, ElemT _default = ElemT());

        // return a reference of the cell (row, col)
        ElemT& get(int row, int col);

        // return the value of the cell (row, col)
        ElemT const_get(int row, int col) const;
    };
}

namespace HPG {
    template<typename ElemT>
    Table<ElemT>::Table(int _rows, int _cols, ElemT _default) : rows(_rows), cols(_cols) {
        table_linearisation.resize(rows * cols, _default);
    }

    template<typename ElemT>
    inline ElemT& Table<ElemT>::get(int row, int col) {
        return table_linearisation[row * cols + col];
    }

    template<typename ElemT>
    inline ElemT Table<ElemT>::const_get(int row, int col) const {
        return table_linearisation[row * cols + col];
    }
}

#endif