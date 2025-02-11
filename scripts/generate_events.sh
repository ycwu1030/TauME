#! /bin/bash

# Define the process
decay_mode=$1
alpha=$2
output_dir=$3
mkdir -p ${output_dir}
multi=10

MG5_dir="/Users/ycwu/Workingspace/MC-Generators/MG5_aMC_v3_5_4"
delphes_dir="/Users/ycwu/Workingspace/Misc/delphes/build/readers"
WORK_dir="/tmp/MG5PROC"
DATA_TMP_dir="/tmp/MG5DATA"
mkdir -p ${WORK_dir}
mkdir -p ${DATA_TMP_dir}
current_dir=$(pwd)
script_dir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd)

s_alpha=$(echo $alpha | sed 's/\./x/g' | sed 's/-/m/g')
run_tag=${decay_mode}_${s_alpha}_$(date '+%y%m%d%H%M%S')
process=Htautau_${run_tag}
process_pi="e- e+ > z h, ( h > ta+ ta-, ta+ > pi+ vt~, ta- > pi- vt )"
process_rho="e- e+ > z h, ( h > ta+ ta-, ta+ > pi+ pi0 vt~, ta- > pi- pi0 vt )"
process_cmd=""
case $decay_mode in
    pi)
        process_cmd=${process_pi}
        ;;
    rho)
        process_cmd=${process_rho}
        ;;
    *)
        echo "Invalid decay mode"
        exit 1
        ;;
esac
process_dir=${WORK_dir}/$process

# Set event number
eventNum=100000
halfEcm=125


cd ${current_dir}
# Generate the process
rm -rf $process.cmd
echo "import model SM_Higgs_tau_CPV_UFO__taudecay_UFO"            >> $process.cmd
echo "generate ${process_cmd}"                                    >> $process.cmd
echo "output ${process_dir}"                                      >> $process.cmd
echo "exit"                                                       >> $process.cmd
cat $process.cmd                                                  >> $process.log
${MG5_dir}/bin/mg5_aMC $process.cmd                               >> $process.log
echo "Process generated"


# Scan the c.m.energy
for nrun in $(seq 1 1 $multi)
do
	rm -rf $process.cmd
	seed=$(awk 'BEGIN{srand();print int(rand()*1000000)}')
	echo "$alpha  $nrun  $seed"
	runname="run_$nrun"
	echo "generate_events $runname"                               >> $process.cmd
	echo "analysis=OFF"                                           >> $process.cmd
	echo "0"                                                      >> $process.cmd
	echo "set ktau 1"                                             >> $process.cmd
	echo "set alp $alpha"                                         >> $process.cmd
	echo "set iseed $seed"                                        >> $process.cmd
	echo "set run_card nevents $eventNum"                         >> $process.cmd
	echo "set run_card lpp1 0"                                    >> $process.cmd
	echo "set run_card lpp2 0"                                    >> $process.cmd
	echo "set run_card ebeam1 $halfEcm"                           >> $process.cmd
	echo "set run_card ebeam2 $halfEcm"                           >> $process.cmd
	echo "set no_parton_cut"                                      >> $process.cmd
	echo "0"                                                      >> $process.cmd
	${process_dir}/bin/madevent $process.cmd >> $process.results
	cp ${process_dir}/Events/$runname/unweighted_events.lhe.gz ${DATA_TMP_dir}/parton_events_${run_tag}_${nrun}.lhe.gz
    gunzip -k ${DATA_TMP_dir}/parton_events_${run_tag}_${nrun}.lhe.gz
	rm -rf $process.results
	rm -rf ${process_dir}/Events/$runname
	rm -rf ${process_dir}/HTML/$runname
	rm -rf ${process_dir}/HTML/results.pkl
	rm -rf ${process_dir}/crossx.html
done

lhefiles=$(ls ${DATA_TMP_dir}/parton_events_${run_tag}_*.lhe)
${delphes_dir}/DelphesLHEF ${script_dir}/gen_card.tcl ${output_dir}/delphes_events_${run_tag}.root $lhefiles

rm -rf $process.cmd $process.log
rm -rf ${process_dir}
rm -rf ${DATA_TMP_dir}/parton_events_${run_tag}_*

cd ${current_dir}
