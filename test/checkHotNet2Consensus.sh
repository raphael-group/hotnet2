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

num_network_permutations=20
num_heat_permutations=1000
num_cores=20

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

# Create networks and vertex weights.
echo Creating networks and vertex weights...
python $hotnet2/test/createSimulatedData.py \
    -n $num_vertices \
    -m $min_num_edges \
    -i 5 5 5 5 5 10 10 10 10 10 \
    -s $seed \
    -nsn $num_vertices_network_1 $num_vertices_network_2 $num_vertices_network_3 \
    -nsw $num_vertices_weights_1 $num_vertices_weights_2 \
    -ivf $data/networks/1/index_vertex.txt $data/networks/2/index_vertex.txt $data/networks/3/index_vertex.txt \
    -elf $data/networks/1/edge_list.txt    $data/networks/2/edge_list.txt    $data/networks/3/edge_list.txt \
    -vwf $data/weights/1/vertex_weight_1.txt $data/weights/2/vertex_weight_2.txt

# Create influence matrices.
for i in 1 2 3
do
    python $hotnet2/makeNetworkFiles.py \
        -nn network_$i \
        -e $data/networks/$i/edge_list.txt \
        -i $data/networks/$i/index_vertex.txt \
        -o $data/networks/$i \
        -np $num_network_permutations \
        -p network_$i \
        -b $beta \
        -c $num_cores
done

# Run HotNet2 consensus
python $hotnet2/HotNet2.py \
    -nf $data/networks/1/network_1_ppr_$beta.h5 \
        $data/networks/2/network_2_ppr_$beta.h5 \
        $data/networks/3/network_3_ppr_$beta.h5 \
    -pnp $data/networks/1/permuted/network_1_ppr_${beta}_##NUM##.h5 \
        $data/networks/2/permuted/network_2_ppr_${beta}_##NUM##.h5 \
        $data/networks/3/permuted/network_3_ppr_${beta}_##NUM##.h5 \
    -hf $data/weights/1/vertex_weight_1.txt \
        $data/weights/2/vertex_weight_2.txt \
    -np $num_network_permutations \
    -hp $num_heat_permutations \
    -o $results/ \
    -c $num_cores \
    --verbose 1

# Compare results with previously computed results.
python $hotnet2/test/compare_consensus_results.py \
    -p $results/consensus/subnetworks.tsv \
    -r referenceConsensusResults.txt
