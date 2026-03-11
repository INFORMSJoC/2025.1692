// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains all the functions that are useful for parsing
    the regex-like string describing the structure of the cutting stages.
    The result is an automaton-like graph.
*/

#ifndef STAGE_GRAMMAR_H
#define STAGE_GRAMMAR_H

#include "../common/hpg.h"

namespace HPG {
    struct StageGrammarState {
        // in this state, we can cut orthogonally to the axes in this container
        std::vector<int> axes;

        // if true, this stage has a single axis and it is the first cutting stage of this axis 
        bool first_stage_optim = false;
        
        // if true, this stage has a single axis and it is the last cutting stage of this axis
        bool last_stage_optim = false;
    };

    struct StageGrammarTransition {
        int input, output;
    };

    struct StageGrammar {
        // for each state of the automaton-like graph 
        std::vector<StageGrammarState> states;

        // transitions in the automaton-like graph
        std::vector<StageGrammarTransition> transitions;

        // state of the initial plate
        int start;

        // state where we can sell items
        int end;

    private:
        // recursively parse an entire expression and returns the input and output of the corresponding gadget. Example : 01(0+2)+10[02]
        StageGrammarTransition parse_expr(std::string& str, size_t& cur_pos);

        // recursively parse a sequence of states and returns the input and output of the corresponding gadget. Example : 01[02]1
        StageGrammarTransition parse_seq(std::string& str, size_t& cur_pos);

        // parse a single state and returns the input and output of the corresponding gadget. Example : [01]
        StageGrammarTransition parse_state(std::string& str, size_t& cur_pos);

    public:
        // default constructor
        StageGrammar();

        // create the automaton-like graph from a regex-like string
        StageGrammar(std::string str);
    };
}

namespace HPG {
    StageGrammar::StageGrammar() {}

    StageGrammar::StageGrammar(std::string str) {
        size_t cur_pos = 0;
        StageGrammarTransition chain = parse_expr(str, cur_pos);
        start = chain.input; 
        end = chain.output;

        // compute for each state if we can apply first or last stage optimization
        for(int state_id = 0;state_id < states.size();state_id++) {
            if(states[state_id].axes.size() != 1) continue;
            int axis = states[state_id].axes[0];
            
            std::vector<bool> is_ancestor(states.size(), false), is_descendent(states.size(), false);
            is_ancestor[state_id] = true;
            is_descendent[state_id] = true;

            for(int _ = 0;_ < states.size();_++) {
                for(const auto& transition : transitions) {
                    if(is_ancestor[transition.output]) is_ancestor[transition.input] = true;
                    if(is_descendent[transition.input]) is_descendent[transition.output] = true;
                }
            }

            states[state_id].first_stage_optim = true;
            states[state_id].last_stage_optim = true;
            for(int other_id = 0;other_id < states.size();other_id++) {
                if(other_id == state_id) continue;

                if(std::count(states[other_id].axes.begin(), states[other_id].axes.end(), axis)) {
                    if(is_ancestor[other_id]) states[state_id].first_stage_optim = false;
                    if(is_descendent[other_id]) states[state_id].last_stage_optim = false;
                }
            }
        }
    }

    StageGrammarTransition StageGrammar::parse_expr(std::string& str, size_t& cur_pos) {
        std::vector<StageGrammarTransition> parallel_gadgets;
        while(true) {
            parallel_gadgets.push_back(parse_seq(str, cur_pos));

            if(cur_pos >= str.size() || str[cur_pos] == ')') break;

            if(str[cur_pos] != '+') {
                std::cerr << "stage grammar error" << std::endl;
                exit(0);
            }
            cur_pos++;
        }

        if(parallel_gadgets.size() == 1) return parallel_gadgets[0];

        int start_id = states.size(); states.push_back({});
        int end_id = states.size(); states.push_back({});
        for(const auto& gadget : parallel_gadgets) {
            transitions.push_back({start_id, gadget.input});
            transitions.push_back({gadget.output, end_id});
        }

        return {start_id, end_id};
    }

    StageGrammarTransition StageGrammar::parse_seq(std::string& str, size_t& cur_pos) {
        StageGrammarTransition head;
        if(str[cur_pos] == '[') {
            cur_pos++; head = parse_state(str, cur_pos);
        }
        else if(str[cur_pos] == '(') {
            cur_pos++; head = parse_expr(str, cur_pos);
            if(str[cur_pos] != ')') {
                std::cerr << "stage grammar error" << std::endl;
                exit(0);
            }
            cur_pos++;
        }
        else if(isdigit(str[cur_pos])) {
            StageGrammarState state;
            state.axes.push_back((int)(str[cur_pos] - '0'));
            int state_id = states.size();
            states.push_back(state);
            cur_pos++;
            head = {state_id, state_id};
        } else {
            std::cerr << "stage grammar error" << std::endl;
            exit(0);
        }

        if(cur_pos >= str.size() || str[cur_pos] == ')' || str[cur_pos] == '+') return head;

        StageGrammarTransition tail = parse_seq(str, cur_pos);
        transitions.push_back({head.output, tail.input});
        return {head.input, tail.output};
    }

    StageGrammarTransition StageGrammar::parse_state(std::string& str, size_t& cur_pos) {
        StageGrammarState state;
        while(isdigit(str[cur_pos])) {
            state.axes.push_back((int)(str[cur_pos] - '0'));
            cur_pos++;
        }
        int state_id = states.size();
        states.push_back(state);

        if(str[cur_pos] != ']') {
            std::cerr << "stage grammar error" << std::endl;
            exit(0);
        }
        cur_pos++;
        return {state_id, state_id};
    }
}

#endif