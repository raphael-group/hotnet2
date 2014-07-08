Consensus Finding
==========================

The file identify_consensus.py can be used to identify consensus subnetworks for general results from different runs of hotnet2. As a result it was slightly modified from mdml's original implementation. The significant change is that it treats each results file as an independent network. So in the case of comparing mutsig and frequency scores, this would cause the same network (ie iRef) to contribute 2 to the edge weight in the case that an edge was reported in both mutsig and frequency. The arguments are also slightly changed.


Arguments
----------------------------
--result_directories/-n 		the location of the actual results directories down to specific deltas ie: hn2example/results/iref/delta_1/

--output_file/-o     output file
-ms/--min_cc_size    minimum connected component size. Default is 2