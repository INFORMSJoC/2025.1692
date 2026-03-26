#/bin/bash

benchmark_name="$1"
instances="$2"
shift 2
options="$@"

rm -rf $benchmark_name
mkdir $benchmark_name
echo $options > $benchmark_name/$benchmark_name.log

ls $instances | sort | parallel -j4 '
    ulimit -v 8000000
    base=`basename {}`
    echo "===== $base ====="
    timeout 7200s ./gkp_gurobi.exe '"$options"' {} '"$benchmark_name"'/$base.hpg 2> '"$benchmark_name"'/$base.log
    ./viz.exe '"$benchmark_name"'/$base.hpg '"$benchmark_name"'/$base.html
    rm -rf '"$benchmark_name"'/$base.hpg
'
