#!/bin/bash

decay_mode=$1
# dflist=(0.0 100.0 -100.0 200.0 -200.0 500.0 -500.0)
dflist=(0.0 200.0 -200.0 500.0 -500.0)
datafolder="/Users/ycwu/Workingspace/WorkingData/TauOOData/STCF_EDM"

script_dir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

for df in ${dflist[@]}
do
    ${script_dir}/generate_events_EDM.sh $decay_mode $df ${datafolder}
done
