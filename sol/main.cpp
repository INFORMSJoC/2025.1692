// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

#include "exact.h"

int main(int argc, char* argv[]) {
    std::ios_base::sync_with_stdio(false);

    if(argc < 3) {
        std::cout << "usage: ./sol.exe input.hpg output.hpg" << std::endl;
        exit(0);
    }

    HPG::HPG hypergraph(argv[argc - 2]);   
    auto sol = HPG::solve(hypergraph);
    save_hpg(sol, argv[argc - 1]);
    return 0;
}