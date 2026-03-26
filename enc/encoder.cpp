// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

#include "encoder_arguments.h"
#include "geometric_instance_file.h"

int main(int argc, char* argv[]) {
    auto start_time = std::chrono::high_resolution_clock::now();

    // parse the command line arguments
    HPG::EncoderArguments arguments(argc, argv);

    // read the geometric instance file
    HPG::GeometricInstanceFile geometric_instance_file(arguments);

    // create a decision hypergraph from the geometric instance
    geometric_instance_file.hpg();

    // write the hypergraph to the output
    save_hpg(geometric_instance_file, arguments.output);

    // compute and print statistics
    size_t nb_vertices = geometric_instance_file.nb_vertices();
    size_t nb_hyperarcs = HPG::nb_hyperarcs(geometric_instance_file);

    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);

    std::cerr << "hypergraph with " << nb_vertices << " vertices and " << nb_hyperarcs << " hyperarcs generated in " << total_time.count() << std::endl;
    return 0;
}