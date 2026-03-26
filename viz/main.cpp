// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

#include "scene.h"

int main(int argc, char* argv[]) {
    if(argc != 3) {
        std::cout << "usage: ./viz.exe input.hpg output.html" << std::endl;
        exit(0);
    }

    HPG::HPG hypergraph(argv[1]);
    HPG::Scene scene(hypergraph);

    std::ifstream viz_template("visualizer_template.html");
    if(!viz_template.is_open()) {
        std::cout << "visualizer_template.html not found" << std::endl;
        exit(0);
    }

    std::ofstream html_output(argv[2]);
    if(!html_output.is_open()) {
        std::cout << "cannot write to " << argv[2] << std::endl;
        exit(0);
    }

    std::string line;
    while(std::getline(viz_template, line)) {
        if(line == "[[DATA]]") {
            scene.write_to(html_output);
        } else {
            html_output << line << std::endl;
        }
    }

    viz_template.close();
    html_output.close();
    return 0;
}