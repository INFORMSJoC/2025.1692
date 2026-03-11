// SPDX-License-Identifier: MPL-2.0
// Copyright (c) 2025 Université de Bordeaux, France
//
// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.

/*
    This file contains the interactions with the general-purpose LP solver.
    This version can use Gurobi, Soplex or Highs.
    You can edit this file if you want to support another solver.
*/

#ifndef MODEL_H
#define MODEL_H

#include "lp_solution.h"
#include "cut_generator.h"

#ifdef USE_GUROBI
#include "gurobi_c++.h"

namespace HPG {
    class Model {
        // Gurobi environment
        GRBEnv* env;

        // Gurobi model
        GRBModel* model;

        // for each item, its variable : the number of times it is produced
        std::vector<GRBVar> item_variable;

        // for each vertex, its flow conservation constraint
        std::map<VertexID, GRBConstr> vertex_constraint;

        // for each hyperarc, its variable : the number of times it is used
        std::map<Hyperarc, GRBVar> hyperarc_variable;

        // for each cut, the corresponding constraint in the model
        std::vector<GRBConstr> cut_constraint;

        // create the vertex if it does not exist, such that the excess flow is at most max_excess
        // and return the corresponding constraint.
        GRBConstr get_vertex(VertexID vertex, float max_excess = 0);

    public:
        // for each vertex u and for each item i, the maximum number of times the vertex u can produce the item i.
        Table<Cnt> maximum_production;

        // for each cut and each vertex, true if the vertex and the item are on the same side
        std::vector<std::vector<bool>> cut_sets;

        // for each cut, its item
        std::vector<int> cut_items;

        // initialize an empty model from a hypergraph instance
        template<typename InstanceT>
        Model(const InstanceT& instance, Table<Cnt>& _maximum_production);

        // delete the model
        ~Model();

        // reoptimize the model using :
        // - the primal simplex algorithm if primal is set to true
        // - the dual simplex algorithm if primal is set to false
        void optimize(bool primal = true);

        // create a new hyperarc
        void new_hyperarc(Hyperarc hyperarc);

        // create a new cut
        void new_cut(int cut_item, const std::vector<bool>& cut_set);

        // return the total profit
        double obj() const;

        // return the number of hyperarcs in the model
        size_t nb_hyperarcs() const;

        // return the dual value of the item's bounds
        double get_item_dual(int item_id) const;

        // return the dual value of the cut
        double get_cut_dual(int cut_id) const;

        // return the optimal fractional solution
        LpSolution get_solution();
    };
}

namespace HPG {
    GRBConstr Model::get_vertex(VertexID vertex, float max_excess) {
        auto it_constr = vertex_constraint.find(vertex);
        if(it_constr != vertex_constraint.end()) return it_constr->second;
        return vertex_constraint[vertex] = model->addConstr(GRBLinExpr(0), GRB_LESS_EQUAL, GRBLinExpr(max_excess));
    }

    template<typename InstanceT>
    Model::Model(const InstanceT& instance, Table<Cnt>& _maximum_production) 
     : maximum_production(_maximum_production) {
        // initialize Gurobi's environment
        env = new GRBEnv(true);
        env->set(GRB_IntParam_LogToConsole, 0);
        env->start();

        // initialize Gurobi's model
        model = new GRBModel(env);

        // create the root vertex
        get_vertex(instance.root, 1);

        // create the items
        for(size_t item_id = 0;item_id < instance.nb_items();item_id++) {
            GRBColumn col;
            col.addTerm(1, get_vertex(instance.items[item_id].vrt));
            item_variable.push_back(
                model->addVar(0, instance.items[item_id].ub, -instance.items[item_id].obj, GRB_CONTINUOUS, col)
            );
        }
    }

    Model::~Model() {
        delete model;
        delete env;
    }

    void Model::optimize(bool primal) {
        if(primal) model->set(GRB_IntParam_Method, 0);
        else model->set(GRB_IntParam_Method, 1);
        model->optimize();
    }

    void Model::new_hyperarc(Hyperarc hyperarc) {
        // check if the hyperarc already exists
        if(hyperarc_variable.find(hyperarc) == hyperarc_variable.end()) {
            // create the new variable
            GRBColumn col;

            // we include the variable in the flow conservation constraints
            col.addTerm(1, get_vertex(hyperarc.r));
            if(hyperarc.u == hyperarc.v) {
                col.addTerm(-2, get_vertex(hyperarc.u));
            } else {
                col.addTerm(-1, get_vertex(hyperarc.u));
                col.addTerm(-1, get_vertex(hyperarc.v));
            }

            // we include the variable in the cuts
            for(int cut_id = 0;cut_id < (int)cut_sets.size();cut_id++) {
                if(!cut_sets[cut_id][hyperarc.r] && (cut_sets[cut_id][hyperarc.u] || cut_sets[cut_id][hyperarc.v])) {
                    col.addTerm(
                        -cut_coeff(
                            cut_items[cut_id], 
                            hyperarc, 
                            cut_sets[cut_id][hyperarc.u], 
                            cut_sets[cut_id][hyperarc.v], 
                            maximum_production
                        ), 
                        cut_constraint[cut_id]
                    );
                }
            }

            // add the variable to the model
            hyperarc_variable[hyperarc] = model->addVar(0, INFINITY, 0, GRB_CONTINUOUS, col);
        }
    }

    void Model::new_cut(int cut_item, const std::vector<bool>& cut_set) {
        GRBLinExpr sum;
        for(auto hyperarc_iterator : hyperarc_variable) {
            Hyperarc hyperarc = hyperarc_iterator.first;
            if(!cut_set[hyperarc.r] && (cut_set[hyperarc.u] || cut_set[hyperarc.v])) {
                sum += cut_coeff(
                    cut_item, 
                    hyperarc, 
                    cut_set[hyperarc.u], 
                    cut_set[hyperarc.v], 
                    maximum_production
                ) * hyperarc_iterator.second;
            }
        }

        cut_constraint.push_back(model->addConstr(item_variable[cut_item] <= sum));
        cut_sets.push_back(cut_set);
        cut_items.push_back(cut_item);
    }

    inline double Model::obj() const {
        return -model->get(GRB_DoubleAttr_ObjVal);
    }

    inline size_t Model::nb_hyperarcs() const {
        return hyperarc_variable.size();
    }

    inline double Model::get_item_dual(int item_id) const {
        return item_variable[item_id].get(GRB_DoubleAttr_RC);
    }

    inline double Model::get_cut_dual(int cut_id) const {
        return cut_constraint[cut_id].get(GRB_DoubleAttr_Pi);
    }

    LpSolution Model::get_solution() {
        LpSolution solution;

        for(int item_id = 0;item_id < item_variable.size();item_id++) {
            solution.item_flow.push_back(
                item_variable[item_id].get(GRB_DoubleAttr_X)
            );
        }

        for(auto hyperarc_iterator : hyperarc_variable) {
            double flow = hyperarc_iterator.second.get(GRB_DoubleAttr_X);
            if(flow > 1e-9) {
                solution.hyperarc_flow[hyperarc_iterator.first] = flow;
            }
        }

        return solution;
    }
}
#endif

#ifdef USE_SOPLEX
#include "soplex.h"

namespace HPG {
    class Model {
        // Soplex model
        soplex::SoPlex* model;

        // current primals
        soplex::DVector primals;

        // current item duals
        soplex::DVector item_duals;

        // current cut duals
        soplex::DVector cut_duals;

        // for each item, its variable : the number of times it is produced
        std::vector<int> item_variable;

        // for each vertex, its flow conservation constraint
        std::map<VertexID, int> vertex_constraint;

        // for each hyperarc, its variable : the number of times it is used
        std::map<Hyperarc, int> hyperarc_variable;

        // for each cut, the corresponding constraint in the model
        std::vector<int> cut_constraint;

        // create the vertex if it does not exist, such that the excess flow is at most max_excess
        // and return the corresponding constraint.
        int get_vertex(VertexID vertex, float max_excess = 0);

    public:
        // for each vertex u and for each item i, the maximum number of times the vertex u can produce the item i.
        Table<Cnt> maximum_production;

        // for each cut and each vertex, true if the vertex and the item are on the same side
        std::vector<std::vector<bool>> cut_sets;

        // for each cut, its item
        std::vector<int> cut_items;

        // initialize an empty model from a hypergraph instance
        template<typename InstanceT>
        Model(const InstanceT& instance, Table<Cnt>& _maximum_production);

        // delete the model
        ~Model();

        // reoptimize the model using :
        // - the primal simplex algorithm if primal is set to true
        // - the dual simplex algorithm if primal is set to false
        void optimize(bool primal = true);

        // create a new hyperarc
        void new_hyperarc(Hyperarc hyperarc);

        // create a new cut
        void new_cut(int cut_item, const std::vector<bool>& cut_set);

        // return the total profit
        double obj() const;

        // return the number of hyperarcs in the model
        size_t nb_hyperarcs() const;

        // return the dual value of the item's bounds
        double get_item_dual(int item_id) const;

        // return the dual value of the cut
        double get_cut_dual(int cut_id) const;

        // return the optimal fractional solution
        LpSolution get_solution();
    };
}

namespace HPG {
    int Model::get_vertex(VertexID vertex, float max_excess) {
        auto it_constr = vertex_constraint.find(vertex);
        if(it_constr != vertex_constraint.end()) return it_constr->second;
        soplex::DSVector row;
        model->addRowReal(soplex::LPRow(-INFINITY, row, max_excess));
        return vertex_constraint[vertex] = model->numRows() - 1;
    }

    template<typename InstanceT>
    Model::Model(const InstanceT& instance, Table<Cnt>& _maximum_production) 
     : maximum_production(_maximum_production) {
        // initialize Soplex's model
        model = new soplex::SoPlex();
        model->setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MINIMIZE);
        model->setIntParam(soplex::SoPlex::VERBOSITY, 0);

        // create the root vertex
        get_vertex(instance.root, 1);

        // create the items
        for(size_t item_id = 0;item_id < instance.nb_items();item_id++) {
            soplex::DSVector col;
            col.add(get_vertex(instance.items[item_id].vrt), 1);
            model->addColReal(soplex::LPCol(-instance.items[item_id].obj, col, instance.items[item_id].ub, 0));
            item_variable.push_back(model->numCols() - 1);
        }
    }

    Model::~Model() {
        delete model;
    }

    void Model::optimize(bool primal) {
        if(primal) model->setIntParam(soplex::SoPlex::ALGORITHM, soplex::SoPlex::ALGORITHM_PRIMAL);
        else model->setIntParam(soplex::SoPlex::ALGORITHM, soplex::SoPlex::ALGORITHM_DUAL);
        
        model->optimize();

        primals.reSize(model->numCols());
        model->getPrimalReal(&primals[0], model->numCols());

        item_duals.reSize(model->numCols());
        model->getRedCostReal(&item_duals[0], model->numCols());

        cut_duals.reSize(model->numRows());
        model->getDualReal(&cut_duals[0], model->numRows());
    }

    void Model::new_hyperarc(Hyperarc hyperarc) {
        // check if the hyperarc already exists
        if(hyperarc_variable.find(hyperarc) == hyperarc_variable.end()) {
            // create the new variable
            soplex::DSVector col;

            // we include the variable in the flow conservation constraints
            col.add(get_vertex(hyperarc.r), 1);
            if(hyperarc.u == hyperarc.v) {
                col.add(get_vertex(hyperarc.u), -2);
            } else {
                col.add(get_vertex(hyperarc.u), -1);
                col.add(get_vertex(hyperarc.v), -1);
            }

            // we include the variable in the cuts
            for(int cut_id = 0;cut_id < (int)cut_sets.size();cut_id++) {
                if(!cut_sets[cut_id][hyperarc.r] && (cut_sets[cut_id][hyperarc.u] || cut_sets[cut_id][hyperarc.v])) {
                    col.add(
                        cut_constraint[cut_id],
                        cut_coeff(
                            cut_items[cut_id], 
                            hyperarc, 
                            cut_sets[cut_id][hyperarc.u], 
                            cut_sets[cut_id][hyperarc.v], 
                            maximum_production
                        )
                    );
                }
            }

            // add the variable to the model
            model->addColReal(soplex::LPCol(0, col, INFINITY, 0));
            hyperarc_variable[hyperarc] = model->numCols() - 1;
        }
    }

    void Model::new_cut(int cut_item, const std::vector<bool>& cut_set) {
        soplex::DSVector constraint;

        for(auto hyperarc_iterator : hyperarc_variable) {
            Hyperarc hyperarc = hyperarc_iterator.first;
            if(!cut_set[hyperarc.r] && (cut_set[hyperarc.u] || cut_set[hyperarc.v])) {
                constraint.add(
                    hyperarc_iterator.second,
                    cut_coeff(
                        cut_item, 
                        hyperarc, 
                        cut_set[hyperarc.u], 
                        cut_set[hyperarc.v], 
                        maximum_production
                    )
                );
            }
        }

        constraint.add(item_variable[cut_item], -1);
        model->addRowReal(soplex::LPRow(0.0, constraint, INFINITY));

        cut_constraint.push_back(model->numRows() - 1);
        cut_sets.push_back(cut_set);
        cut_items.push_back(cut_item);
    }

    inline double Model::obj() const {
        return -model->objValueReal();
    }

    inline size_t Model::nb_hyperarcs() const {
        return hyperarc_variable.size();
    }

    inline double Model::get_item_dual(int item_id) const {
        return item_duals[item_variable[item_id]];
    }

    inline double Model::get_cut_dual(int cut_id) const {
        return -cut_duals[cut_constraint[cut_id]];
    }

    LpSolution Model::get_solution() {
        LpSolution solution;

        for(int item_id = 0;item_id < item_variable.size();item_id++) {
            solution.item_flow.push_back(
                primals[item_variable[item_id]]
            );
        }

        for(auto hyperarc_iterator : hyperarc_variable) {
            double flow = primals[hyperarc_iterator.second];
            if(flow > 1e-9) {
                solution.hyperarc_flow[hyperarc_iterator.first] = flow;
            }
        }

        return solution;
    }
}
#endif

#ifdef USE_HIGHS
#include "Highs.h"

namespace HPG {
    class Model {
        // Highs model
        Highs* model;

        // current solution
        HighsSolution current_solution;

        // for each item, its variable : the number of times it is produced
        std::vector<int> item_variable;

        // for each vertex, its flow conservation constraint
        std::map<VertexID, int> vertex_constraint;

        // for each hyperarc, its variable : the number of times it is used
        std::map<Hyperarc, int> hyperarc_variable;

        // for each cut, the corresponding constraint in the model
        std::vector<int> cut_constraint;

        // create the vertex if it does not exist, such that the excess flow is at most max_excess
        // and return the corresponding constraint.
        int get_vertex(VertexID vertex, float max_excess = 0);

    public:
        // for each vertex u and for each item i, the maximum number of times the vertex u can produce the item i.
        Table<Cnt> maximum_production;

        // for each cut and each vertex, true if the vertex and the item are on the same side
        std::vector<std::vector<bool>> cut_sets;

        // for each cut, its item
        std::vector<int> cut_items;

        // initialize an empty model from a hypergraph instance
        template<typename InstanceT>
        Model(const InstanceT& instance, Table<Cnt>& _maximum_production);

        // delete the model
        ~Model();

        // reoptimize the model using :
        // - the primal simplex algorithm if primal is set to true
        // - the dual simplex algorithm if primal is set to false
        void optimize(bool primal = true);

        // create a new hyperarc
        void new_hyperarc(Hyperarc hyperarc);

        // create a new cut
        void new_cut(int cut_item, const std::vector<bool>& cut_set);

        // return the total profit
        double obj() const;

        // return the number of hyperarcs in the model
        size_t nb_hyperarcs() const;

        // return the dual value of the item's bounds
        double get_item_dual(int item_id) const;

        // return the dual value of the cut
        double get_cut_dual(int cut_id) const;

        // return the optimal fractional solution
        LpSolution get_solution();
    };
}

namespace HPG {
    int Model::get_vertex(VertexID vertex, float max_excess) {
        auto it_constr = vertex_constraint.find(vertex);
        if(it_constr != vertex_constraint.end()) return it_constr->second;
        model->addRow(-INFINITY, max_excess, 0, nullptr, nullptr);
        return vertex_constraint[vertex] = model->getNumRow() - 1;
    }

    template<typename InstanceT>
    Model::Model(const InstanceT& instance, Table<Cnt>& _maximum_production) 
     : maximum_production(_maximum_production) {
        // initialize Highs's model
        model = new Highs();
        model->changeObjectiveSense(ObjSense::kMinimize);
        model->setOptionValue("log_to_console", false);

        // create the root vertex
        get_vertex(instance.root, 1);

        // create the items
        for(size_t item_id = 0;item_id < instance.nb_items();item_id++) {
            std::vector<int> indices = {get_vertex(instance.items[item_id].vrt)}; 
            std::vector<double> coeffs = {1};
            model->addCol(-instance.items[item_id].obj, 0, instance.items[item_id].ub, 1, &indices[0], &coeffs[0]);
            item_variable.push_back(model->getNumCol() - 1);
        }
    }

    Model::~Model() {
        delete model;
    }

    void Model::optimize(bool primal) {
        if(primal) {
            model->setOptionValue("solver", "simplex");
            model->setOptionValue("simplex_strategy", 4);
        }
        else {
            model->setOptionValue("solver", "simplex");
            model->setOptionValue("simplex_strategy", 1);
        }

        model->run();
        current_solution = model->getSolution();
    }

    void Model::new_hyperarc(Hyperarc hyperarc) {
        // check if the hyperarc already exists
        if(hyperarc_variable.find(hyperarc) == hyperarc_variable.end()) {
            // create the new variable
            std::vector<int> indices; std::vector<double> coeffs;

            // we include the variable in the flow conservation constraints
            indices.push_back(get_vertex(hyperarc.r)); coeffs.push_back(1);
            if(hyperarc.u == hyperarc.v) {
                indices.push_back(get_vertex(hyperarc.u)); coeffs.push_back(-2);
            } else {
                indices.push_back(get_vertex(hyperarc.u)); coeffs.push_back(-1);
                indices.push_back(get_vertex(hyperarc.v)); coeffs.push_back(-1);
            }

            // we include the variable in the cuts
            for(int cut_id = 0;cut_id < (int)cut_sets.size();cut_id++) {
                if(!cut_sets[cut_id][hyperarc.r] && (cut_sets[cut_id][hyperarc.u] || cut_sets[cut_id][hyperarc.v])) {
                    indices.push_back(cut_constraint[cut_id]);
                    coeffs.push_back(
                        cut_coeff(
                            cut_items[cut_id], 
                            hyperarc, 
                            cut_sets[cut_id][hyperarc.u], 
                            cut_sets[cut_id][hyperarc.v], 
                            maximum_production
                        )
                    );
                }
            }

            // add the variable to the model
            model->addCol(0, 0, INFINITY, indices.size(), &indices[0], &coeffs[0]);
            hyperarc_variable[hyperarc] = model->getNumCol() - 1;
        }
    }

    void Model::new_cut(int cut_item, const std::vector<bool>& cut_set) {
        std::vector<int> indices; std::vector<double> coeffs;

        for(auto hyperarc_iterator : hyperarc_variable) {
            Hyperarc hyperarc = hyperarc_iterator.first;
            if(!cut_set[hyperarc.r] && (cut_set[hyperarc.u] || cut_set[hyperarc.v])) {
                indices.push_back(hyperarc_iterator.second);
                coeffs.push_back(
                    cut_coeff(
                        cut_item, 
                        hyperarc, 
                        cut_set[hyperarc.u], 
                        cut_set[hyperarc.v], 
                        maximum_production
                    )
                );
            }
        }

        indices.push_back(item_variable[cut_item]);
        coeffs.push_back(-1);

        model->addRow(0, INFINITY, indices.size(), &indices[0], &coeffs[0]);

        cut_constraint.push_back(model->getNumRow() - 1);
        cut_sets.push_back(cut_set);
        cut_items.push_back(cut_item);
    }

    inline double Model::obj() const {
        return -model->getObjectiveValue();
    }

    inline size_t Model::nb_hyperarcs() const {
        return hyperarc_variable.size();
    }

    inline double Model::get_item_dual(int item_id) const {
        return current_solution.col_dual[item_variable[item_id]];
    }

    inline double Model::get_cut_dual(int cut_id) const {
        return -current_solution.row_dual[cut_constraint[cut_id]];
    }

    LpSolution Model::get_solution() {
        LpSolution solution;

        for(int item_id = 0;item_id < item_variable.size();item_id++) {
            solution.item_flow.push_back(
                current_solution.col_value[item_variable[item_id]]
            );
        }

        for(auto hyperarc_iterator : hyperarc_variable) {
            double flow = current_solution.col_value[hyperarc_iterator.second];
            if(flow > 1e-9) {
                solution.hyperarc_flow[hyperarc_iterator.first] = flow;
            }
        }

        return solution;
    }
}
#endif

#endif