# README #

### DESCRIPTION ###
The `checkHotNet2Consensus.sh` script checks the output of the HotNet2 consensus against previously computed results on simulated data.  To run it with default settings, simply run `sh checkHotNet2Consesus.sh`.  From start to finish, the script finishes in under 8 minutes using less than 1 GB of memory and 250 MB of temporary storage space on a modern eight-core machine with the default settings.

The `createSimulatedData.py` script creates a random graph with implanted cliques with high vertex weights.  The `exampleConsensusResults.txt` file contains the example consensus results under the default settings. The results should be highly similar to `exampleConsensusResults.txt`, but they may not be identical.
