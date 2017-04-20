#!/usr/bin/env bash

# These commands run HotNet2 on the processed network and mutation data used for our Nature Genetics paper.

# Set location of HotNet2, number of cores, number of network permutations, and number of heat permutations.
hotnet2=..
num_cores=-1
num_network_permutations=100
num_heat_permutations=1000

# Create network data.
python $hotnet2/makeNetworkFiles.py \
    -e  data/networks/hint+hi2012/hint+hi2012_edge_list \
    -i  data/networks/hint+hi2012/hint+hi2012_index_gene \
    -nn hint+hi2012 \
    -p  hint+hi2012 \
    -b  0.4 \
    -o  data/networks/hint+hi2012 \
    -np $num_network_permutations \
    -c  $num_cores

python $hotnet2/makeNetworkFiles.py \
    -e  data/networks/irefindex9/irefindex9_edge_list \
    -i  data/networks/irefindex9/irefindex9_index_gene \
    -nn irefindex9 \
    -p  irefindex9 \
    -b  0.45 \
    -o  data/networks/irefindex9 \
    -np $num_network_permutations \
    -c  $num_cores

python $hotnet2/makeNetworkFiles.py \
    -e  data/networks/multinet/multinet_edge_list \
    -i  data/networks/multinet/multinet_index_gene \
    -nn multinet \
    -p  multinet \
    -b  0.5 \
    -o  data/networks/multinet \
    -np $num_network_permutations \
    -c  $num_cores

# Create heat data.
python $hotnet2/makeHeatFile.py \
    scores \
    -hf data/heats/pan12.gene2freq.txt \
    -o  data/heats/pan12.gene2freq.json \
    -n  pan12.freq

python $hotnet2/makeHeatFile.py \
    scores \
    -hf data/heats/pan12.gene2mutsig.txt \
    -o  data/heats/pan12.gene2mutsig.json \
    -n  pan12.mutsig

# Run HotNet2.
python $hotnet2/HotNet2.py \
    -nf  data/networks/hint+hi2012/hint+hi2012_ppr_0.4.h5 \
         data/networks/irefindex9/irefindex9_ppr_0.45.h5 \
         data/networks/multinet/multinet_ppr_0.5.h5 \
    -pnp data/networks/hint+hi2012/permuted/hint+hi2012_ppr_0.4_##NUM##.h5 \
         data/networks/irefindex9/permuted/irefindex9_ppr_0.45_##NUM##.h5 \
         data/networks/multinet/permuted/multinet_ppr_0.5_##NUM##.h5 \
    -hf  data/heats/pan12.gene2freq.json \
         data/heats/pan12.gene2mutsig.json \
    -np  $num_network_permutations \
    -hp  $num_heat_permutations \
    -o   results \
    -c   $num_cores
