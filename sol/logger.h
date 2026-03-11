// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains the functions to log data of the main algorithm.
*/

#ifndef LOGGER_H
#define LOGGER_H

#include "../common/types.h"

namespace HPG {
    struct Logger {
    private:
        // starting time of the 
        std::chrono::time_point<std::chrono::steady_clock> start_time;

        // logger's mutex for multithreading
        std::mutex mtx;

        // print the prefix time [lb, ub] :
        void prefix();

    public:
        // number of hyperarc generation rounds in the current cutting plane round
        int nb_lp_iters = 0;

        // maximum number of labels in a vertex during the current labeling round
        size_t max_size = 0;

        // total number of labels during the current labeling round
        size_t total_size = 0;
        
        // number of hyperarcs with positive flow in the current fractional solution
        size_t active_hyperarcs = 0;
        
        // number of generated hyperarcs
        size_t nb_hyperarcs = 0;
        
        // number of generated cuts 
        size_t nb_cuts = 0;

        // current LP value
        Obj lp_ub;

        // current upper bound
        Obj ub;

        // current lower bound
        Obj lb;

        // verbosity of the logger
        // true : print everything
        // false : print nothing
        bool verbose = true;

        // reset the global parameters
        void start();

        // log an event
        void event(std::string str);

        // log an lp event : the end of a cutting plane round
        void lp_event(bool restart = false);

        // log a label event : the end of a labeling round
        void label_event();
    }; 

    // the global logger
    Logger logger;
}

namespace HPG {
    void Logger::start() {
        start_time = std::chrono::steady_clock::now();
        lb = 0;
        ub = INFINITY;
    }

    void Logger::prefix() {
        auto end_time = std::chrono::steady_clock::now();
        std::cerr << 
            std::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()) 
            << " [" << (int)(lb + 1e-3) << ", " << (int)(ub + 1e-3) << "]: ";
    }

    void Logger::event(std::string str) {
        std::lock_guard _(mtx);
        if(verbose) {
            prefix();
            std::cerr << str << std::endl;
        }
    }

    void Logger::lp_event(bool restart) {
        std::string str = "LP\n\tbnd = " + std::to_string(lp_ub) + 
        "\n\tits = " + std::to_string(nb_lp_iters) + 
        "\n\thyp = " + std::to_string(active_hyperarcs) + " / " + std::to_string(nb_hyperarcs) + 
        "\n\tcut = " + std::to_string(nb_cuts);
        if(restart) str += "\n\trestart";
        event(str);
    }

    void Logger::label_event() {
        event("LABEL\n\ttot = " + std::to_string(total_size) + "\n\tmax = " + std::to_string(max_size));
    }
}

#endif