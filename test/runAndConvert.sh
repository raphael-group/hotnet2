#!/bin/sh

# Perform the runs
rm -r tmp
rm -r ../../hotnet2-test/test/tmp-converted
sh checkHotNet2Consensus.sh

python convertSimulatedDataToOldFormat.py \
	-nf tmp/data/networks/1/original/network_1_ppr_0.4.h5 \
	tmp/data/networks/2/original/network_2_ppr_0.4.h5 \
	tmp/data/networks/3/original/network_3_ppr_0.4.h5 \
	-hf tmp/data/weights/1/vertex_weight_1.txt \
	tmp/data/weights/2/vertex_weight_2.txt \
	-o ../../hotnet2-test/test/tmp-converted/data/

cd ../../hotnet2-test/test
sh checkHotNet2ConsensusConverted.sh

# Check the single runs to be the same
cd ../../hotnet2/test
python compare_subnetworks.py -d tmp/results/network_1-vertex_weight_1/ ../../hotnet2-test/test/tmp-converted/results/network_1_weight_1/
python compare_subnetworks.py -d tmp/results/network_2-vertex_weight_1/ ../../hotnet2-test/test/tmp-converted/results/network_2_weight_1/
python compare_subnetworks.py -d tmp/results/network_3-vertex_weight_1/ ../../hotnet2-test/test/tmp-converted/results/network_3_weight_1/
python compare_subnetworks.py -d tmp/results/network_1-vertex_weight_2/ ../../hotnet2-test/test/tmp-converted/results/network_1_weight_2/
python compare_subnetworks.py -d tmp/results/network_2-vertex_weight_2/ ../../hotnet2-test/test/tmp-converted/results/network_2_weight_2/
python compare_subnetworks.py -d tmp/results/network_3-vertex_weight_2/ ../../hotnet2-test/test/tmp-converted/results/network_3_weight_2/