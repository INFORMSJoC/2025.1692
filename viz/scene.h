// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

#ifndef SCENE_H
#define SCENE_H

#include "../common/hpg.h"

namespace HPG {
    struct Box {
        std::vector<int> lbs, ubs;
        size_t item;
    };

    class Scene {
        std::vector<int> dims;
        std::vector<Box> boxes;
    public:
        Scene(HPG hpg);
        void write_to(std::ostream& out);
    };
}

namespace HPG {
    std::vector<int> parse_name(std::string name) {
        size_t pos = name.find('/');
        name = std::string(name.begin(), name.begin() + std::min(pos, name.size()));
        std::replace(name.begin(), name.end(), ',', ' ');
        std::stringstream ss;
        ss << name;
        std::vector<int> coords;
        while (true) {
            int coord; ss >> coord;
            if (ss.fail()) break;
            else coords.push_back(coord);
        }
        coords.resize(3, 0);
        return coords;
    }

    Scene::Scene(HPG hpg) {
        std::vector<size_t> iEdge(hpg.nb_vertices(), 0);

        auto create_scene_rec = [&] (auto self, VertexID r, std::vector<int> shifts) -> int {
            auto coords_r = parse_name(hpg.vertex_names[r]);

            int type = -1;
            if (iEdge[r] < hpg.hyperarcs[r].size()) {
                const auto& edge = hpg.hyperarcs[r][iEdge[r]];

                auto coords_u = parse_name(hpg.vertex_names[edge.u]);
                auto coords_v = parse_name(hpg.vertex_names[edge.v]);

                type = std::max(type, self(self, edge.u, shifts));

                for (int dim = 0; dim < 3; dim++) {
                    if (coords_r[dim] != 0 && coords_r[dim] == coords_u[dim] + coords_v[dim]) {
                        shifts[dim] += coords_u[dim];
                    }
                }

                type = std::max(type, self(self, edge.v, shifts));

                iEdge[r] += 1;
            }

            if (type != -1 && std::count(coords_r.begin(), coords_r.end(), 0) != 3) {
                auto lbs = shifts;
                auto ubs = shifts;
                ubs[0] += coords_r[0]; ubs[1] += coords_r[1]; ubs[2] += coords_r[2];
                boxes.push_back({lbs, ubs, (size_t)type});
                return -1;
            }

            for (size_t item = 0; item < hpg.nb_items(); item++) {
                if (r == hpg.items[item].vrt) {
                    return item;
                }
            }
            return type;
        };

        dims = parse_name(hpg.vertex_names[hpg.root]);
        create_scene_rec(create_scene_rec, hpg.root, {0, 0, 0});

        for(auto& box : boxes) {
            box.ubs[0] = std::max(box.ubs[0], box.lbs[0] + 1);
            box.ubs[1] = std::max(box.ubs[1], box.lbs[1] + 1);
            box.ubs[2] = std::max(box.ubs[2], box.lbs[2] + 1);
        }

        for(int dim = 0; dim < 3;dim++) {
            dims[dim] = std::max(dims[dim], 1);
        }
    }

    void Scene::write_to(std::ostream& out) {
        out << "let dims = [" << dims[0] << ", " << dims[1] << ", " << dims[2] << "];" << std::endl;

        out << "let boxes = [";
        for(auto& box : boxes) {
            out << "{ lbs: [" << box.lbs[0] << ", " << box.lbs[1] << ", " << box.lbs[2] << "], "
            << "ubs: [" << box.ubs[0] << ", " << box.ubs[1] << ", " << box.ubs[2] << "], "
            << "item: " << box.item << "},";
        }
        out << "];" << std::endl;
    }
}

#endif
