[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# An Efficient Solver for Integral Flows in Decision Hypergraphs with Applications to Cutting Problems

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MPL 2.0 License](LICENSE.md).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
[An Efficient Solver for Integral Flows in Decision Hypergraphs with Applications to Cutting Problems](https://doi.org/10.1287/ijoc.2025.1692) by A. Léonard and F. Clautiaux.

**Important: This code is being developed on an on-going basis at 
https://gitlab.inria.fr/edge/hypergraphsolver. Please go there if you would like to
get a more recent version or would like support**

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2025.1692

https://doi.org/10.1287/ijoc.2025.1692.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{leonard2026,
  author =        {A. Léonard and F. Clautiaux},
  publisher =     {INFORMS Journal on Computing},
  title =         {{An Efficient Solver for Integral Flows in Decision Hypergraphs with Applications to Cutting Problems}},
  year =          {2026},
  doi =           {10.1287/ijoc.2025.1692.cd},
  url =           {https://github.com/INFORMSJoC/2025.1692},
  note =          {Available for download at https://github.com/INFORMSJoC/2025.1692},
}  
```

## Description

HPG is a solver for guillotine cutting problems using decision hypergraphs.

This repository contains four useful programs :
* `enc` : An encoder that transform many variants of cutting probems into the decision hypergraph framework.
* `sol` : The actual solver that computes the lower and upper bounds on solutions of an instance given as a decision hypergraph.
* `viz` : A simple solution visualizer which can be used for 2D and 3D problems.
* `gkp` : The actual solver for geometric problems, which handles the hypergraph implicitly.

The folder `instances` contains instances from the literature, and the folder `results` contains the results on the instances from the literature.

## Installation guide 

The best performances are obtained using Gurobi solver. By default, if Gurobi is installed on your computer, the solver will use it. If Gurobi is not available to you, the most recent version of HigHs is installed.

### Ubuntu / Debian
Step 1 : install the following packages
```sh
sudo apt update
sudo apt install g++ cmake
```

(Optional, for SoPlex):
```sh
sudo apt install libboost-all-dev
```

(Recommended, for Gurobi): Get a Gurobi license and download Gurobi Optimizer from [their website](https://www.gurobi.com/downloads/).

Step 2 : Compile HPG with cmake.
```sh
mkdir build
cd build
cmake ..
make
cd ..
mv build/*.exe .
mv build/*.html .
```

If everyting worked, you should see four executables : `enc.exe`, `gkp_<LP solver name>.exe`, `sol_<LP solver name>.exe` and `viz.exe`.

## Exemple

Suppose that we want to solve the any-stage guillotine 2-dimensional knapsack instance P1_100_200_25_1 (without rotation) from Set9.

Step 1 : We solve the instance with the following command.
```sh
./gkp_highs.exe --rotations=01 instances/Set9/P1_100_200_25_1 solution.hpg
./gkp_soplex.exe --rotations=01 instances/Set9/P1_100_200_25_1 solution.hpg
./gkp_gurobi.exe --rotations=01 instances/Set9/P1_100_200_25_1 solution.hpg
```

Step 2 : We create the HTML solution file using the visualizer.
```sh
./viz.exe solution.hpg solution.html
```

Step 3 : We open the HTML solution with your favorite browser.
```sh
firefox solution.html
```
In the recommended procedure described above, the hypergraph is handled implicitly ; the hyperarcs are never stored in memory.
If you wish to inspect the hypergraph, you can use the encoder with
```sh
./enc.exe --rotations=01 instances/Set9/P1_100_200_25_1 problem.hpg
```

If you want to solve your own hypergraph-formulated problem, you can use 
```sh
./sol_highs.exe problem.hpg solution.hpg
./sol_soplex.exe problem.hpg solution.hpg
./sol_gurobi.exe problem.hpg solution.hpg
```

## Specification of the hpg file format

A hpg file contains 5 sections.
```txt
... ignored comments ...
<symbols>
... list of vertices names ...
<edges>
... specifications of the hyperarcs ...
<root>
... specification of the root ...
<items>
... specifications of the items ...
```

The `<symbols>` section specify a list of names for the vertices. The first name given is the one for vertex ID 0, the second name given is the one for vertex ID 1, and so on. This section is useful for the visualizer software, and to link vertices from the original graph and vertices in the solution.
The `<edges>` section is divided into multiple chunks. Each chunk starts with a single vertex ID `t` on a line.
Every other line in the section contains two Vertex ID `s1` and `s2`, indicating an hyperarc with target `t` and sources `s1` and `s2`.
The `<root>` section should contain a single line containing the vertex ID of the root `r` (ie, the type of the initial plate).
The `<items>` section should contain one line per item containing a vertex ID `i`, a real number `p` and a positive integer `b`,
indicating that plates of type `i` can be sold at most `b` times for a unit profit of `p`.

For example, the example in the paper can be written as follows.
```txt
<symbols>
r
a
b
i0
i1
i2
<edges>
0
1 1
4 2
1
3 4
2
4 5
<root>
0
<items>
3 3 1
4 1 2
5 1 1
```

## Specification of a geometric instance file

Most instances in the literature are geometric instances, where items corresponds to boxes (1D, 2D or 3D).
Suppose that we are in dimension `k`.
The first line of a geometric instance file contains `k` integers, corresponding to the dimensions of the initial plate.
The second line contains a single integer, the number of items.
Follows one line per item, each containing `k + 2` numbers. The first `k` numbers are integers, corresponding to the dimensions
of the item. The `k + 1`-th number is a real number corresponding to the profit of the item, and the last integer corresponds to 
the maximum number of times the item can be sold. 

## Documentation for `enc`

The program `enc` is an encoder that takes as input a geometric instance file, and outputs a corresponding hpg file.
It uses many rules to reduce the generated hypergraphs, using symmetries and subset sums.
The syntax is the following :
```sh
./enc.exe [options] input output
```
The same arguments are used for the direct solver :
```sh
./gkp_gurobi.exe [options] input output
```

The options are the following :
* `-d` or `--dim` : 
Specifies the dimension of the problem. 
Example : `--dim=3`. 
Default : 2.
* `-s` or `--stages` : 
Specifies the cutting stages directions. The cutting stages directions are given as a string in a regex-like grammar. 
Examples : `--stages=2[01]`, `--stages=01+10`, `--stages=[01]2`, ... 
Default : any-stage.
* `-t` or `--trim` : 
Specifies in which dimensions we are allowed to do trimming cuts. 
Example : `--trim=101` allows trimming only in axes 0 and 2.
Default : all trimming directions allowed.
* `-r` or `--rotations` : Specifies the set of rotations that are allowed.
Example : `--rotations=01` forbid the 2D rotations. `--rotations=012,210` only allows rotations around axis 1.
Default : all rotations allowed. 
* `-m` or `--measure` : If specified, forces the profit of an item to become its measure (length, area, volume, ...).
* `-nc` or `--no-cuts` : If specified, disables the valid inequalities.

## Documentation for `viz`

The program `viz` is the visualization program. It takes as input a hpg file describing the solution,
and create an HTML solution file that can be used to visualize in 2D or 3D.
The syntax is the following :
```sh
./viz.exe input output
```
You can use the arrows to rotate the solution.

## License

This project is licensed under the Mozilla Public License 2.0 (MPL-2.0).
