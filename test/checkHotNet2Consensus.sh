#!/usr/bin/env bash

# Setup directories.
hotnet2=../
hotnet2_test=/tmp/hotnet2_consensus_test
data=$hotnet2_test/data
results=$hotnet2_test/results

# Setup parameters.
beta=0.4
num_network_permutations=100
num_heat_permutations=200
num_cores=8

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
    -n 250 \
    -m 750 \
    -i 5 5 5 5 5 10 10 10 10 10 \
    -s 0 \
    -nsn 225 235 245 \
    -nsw 230 240 \
    -ivf $data/networks/1/index_vertex.txt $data/networks/2/index_vertex.txt $data/networks/3/index_vertex.txt \
    -elf $data/networks/1/edge_list.txt    $data/networks/2/edge_list.txt    $data/networks/3/edge_list.txt \
    -vwf $data/weights/1/vertex_weight.txt $data/weights/2/vertex_weight.txt

# Create influence matrices.
echo Creating influence matrices...
for i in 1 2 3
do
    python $hotnet2/bin/createPPRMat.py \
        -e $data/networks/$i/edge_list.txt \
        -i $data/networks/$i/index_vertex.txt \
        -o $data/networks/$i/original \
        -p network_$i \
        -b $beta
done

# Permute networks.
echo Permuting networks...
for i in 1 2 3
do
    python $hotnet2/bin/permuteNetwork.py \
        -e $data/networks/$i/original/network_"$i"_edge_list \
        -o $data/networks/$i/permuted \
        -p network_$i \
        -n $num_network_permutations \
        -c $num_cores
done

# Create permuted influence matrices.
echo Creating permuted influence matrices...
for i in 1 2 3
do
    for k in `seq $num_network_permutations`
    do
        python $hotnet2/bin/createPPRMat.py \
            -e $data/networks/$i/original/network_"$i"_edge_list \
            -i $data/networks/$i/original/network_"$i"_index_genes \
            -o $data/networks/$i/permuted/$k \
            -p network_$i \
            -b $beta
    done
done

# Run individual HotNet2 runs.
echo Running HotNet2...
for i in 1 2 3
do
    for j in 1 2
    do
        python $hotnet2/runHotNet2.py \
            -mf $data/networks/$i/original/network_"$i"_ppr_"$beta".h5 \
            -if $data/networks/$i/original/network_"$i"_index_genes \
            -hf $data/weights/$j/vertex_weight.txt \
            -pnp $data/networks/$i/permuted/##NUM##/network_"$i"_ppr_"$beta".h5 \
            -dp $num_network_permutations \
            -sp $num_heat_permutations \
            -o $results/network_"$i"_weight_"$j" \
            -c $num_cores
    done
done

# Run consensus procedure.
echo Running consensus procedure...
python $hotnet2/identifyConsensus.py \
    -d $results/network_1_weight_1 $results/network_1_weight_2 \
       $results/network_2_weight_1 $results/network_2_weight_2 \
       $results/network_3_weight_1 $results/network_3_weight_2 \
    -n network_1 network_1 network_2 network_2 network_3 network_3 \
    -o $results/consensus_results.txt

# Compare results with previously computed results; currently, compare manually.
diff $results/consensus_results.txt exampleConsensusResults.txt

# Remove data and results.
rm -rf $hotnet2_consensus_test
