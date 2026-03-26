#!/bin/bash

./../benchmark.sh Set1 "../instances/Set1/*"   --rotations=01
./../benchmark.sh Set2 "../instances/Set2/*"   --rotations=01
./../benchmark.sh Set3 "../instances/Set3/*"   --rotations=01
./../benchmark.sh Set4 "../instances/Set4/*"   --rotations=01
./../benchmark.sh Set5 "../instances/Set5/*"   --rotations=01
./../benchmark.sh Set6 "../instances/Set6/*"   --rotations=01
./../benchmark.sh Set7 "../instances/Set7/*"   --rotations=01
./../benchmark.sh Set8 "../instances/Set8/*"   --rotations=01 --no-cuts 
./../benchmark.sh Set9 "../instances/Set9/*"   --rotations=01
./../benchmark.sh Set10 "../instances/Set10/*" --rotations=01
./../benchmark.sh Set11 "../instances/Set11/*" --rotations=01