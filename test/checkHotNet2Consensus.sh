#!/usr/bin/env bash

# Setup directories.
hotnet2=../
hotnet2_test=tmp/
data=$hotnet2_test/data
results=$hotnet2_test/results

# Set parameters.
beta=0.4

num_vertices=250
min_num_edges=750
num_vertices_network_1=225
num_vertices_network_2=235
num_vertices_network_3=245
num_vertices_weights_1=230
num_vertices_weights_2=240
seed=0

num_network_permutations=50
num_heat_permutations=200
num_cores=1

mkdir -p $data
mkdir -p $results

for i in 1 2 3
do
    mkdir -p $data/networks/$i
    mkdir -p $data/networks/$i/permuted
done

for j in 1 2
do
    mkdir -p $data/weights/$j
done

# # Create networks and vertex weights.
# echo Creating networks and vertex weights...
python $hotnet2/test/createSimulatedData.py \
    -n $num_vertices \
    -m $min_num_edges \
    -i 5 5 5 5 5 10 10 10 10 10 \
    -s $seed \
    -nsn $num_vertices_network_1 $num_vertices_network_2 $num_vertices_network_3 \
    -nsw $num_vertices_weights_1 $num_vertices_weights_2 \
    -ivf $data/networks/1/index_vertex.txt $data/networks/2/index_vertex.txt $data/networks/3/index_vertex.txt \
    -elf $data/networks/1/edge_list.txt    $data/networks/2/edge_list.txt    $data/networks/3/edge_list.txt \
    -vwf $data/weights/1/vertex_weight.txt $data/weights/2/vertex_weight.txt

# Create influence matrices.
for i in 1 2 3
do
    python $hotnet2/makeRequiredPPRFiles.py \
        -nn network_$i \
        -e $data/networks/$i/edge_list.txt \
        -i $data/networks/$i/index_vertex.txt \
        -o $data/networks/$i/original \
        -np $num_network_permutations \
        -p network_$i \
        -b $beta \
        -c $num_cores
done

# Run HotNet2 consensus
python $hotnet2/runHotNet2.py \
    -nf $data/networks/1/original/network_1_ppr_$beta.h5 \
        $data/networks/2/original/network_2_ppr_$beta.h5 \
        $data/networks/3/original/network_3_ppr_$beta.h5 \
    -hf $data/weights/{1,2}/vertex_weight.txt \
    -dp $num_network_permutations \
    -sp $num_heat_permutations \
    -o $results/ \
    -c $num_cores

# Compare results with previously computed results; currently, compare manually.
diff computedConsensusResults.txt referenceConsensusResults.txt

# Remove data and results.
# rm -rf $hotnet2_test
