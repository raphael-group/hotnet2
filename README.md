HotNet2
=======================

HotNet is an algorithm for finding significantly altered subnetworks in a large gene interaction
network. While originally developed for use with cancer mutation data, the current release of
HotNet also supports the application in which scores can be assigned to genes in the network
(previously called generalizedHotNet).

Requirements
------------------------

* Linux/Unix
* [Python 2.7](http://python.org/)
* [NumPy 1.6.2](http://www.numpy.org/)
* [SciPy 0.10.1](http://www.scipy.org/)
* [NetworkX 1.7](http://networkx.github.io/)

HotNet will likely work with additional versions of Python, NetworkX, NumPy, and SciPy, but
alternative configurations have not been tested.

Support
------------------------
For support using HotNet, please visit the [HotNet Google Group](https://groups.google.com/forum/#!forum/hotnet-users).

HotNet2 vs. Classic HotNet
------------------------
This distribution contains two related algorithms for finding significantly altered subnetworks in
a large gene interaction network: the original HotNet algorithm "classic HotNet", and an updated
version "HotNet2".  HotNet2 differs from classic HotNet in several important ways.  First, HotNet2 
ses a new heat diffusion kernel analogous to random walk with restart that better captures the
local topology of the interaction network surrounding a protein compared to the general heat
diffusion process used by HotNet.  HotNet2 also uses an asymmetric influence score and different
permutation testing and parameter selection procedures.

For more details, please refer to the publications listed at the end of this README.

Simple runs
------------------------
To get started running HotNet2 quickly and easily, use the `simpleRun.py` Python script.  You must
provide the following parameters:

        =====================================================================================
        | PARAMETER NAME         | DESCRIPTION                                              |
        =====================================================================================
        |-mf/--infmat_file       |Path to .mat file containing influence matrix downloaded  |
        |                        |from http://compbio.cs.brown.edu/projects/hotnet/         |
        -------------------------------------------------------------------------------------
        |-if/--infmat_index_file |Path to gene-index mapping file downloaded from           |
        |                        |http://compbio.cs.brown.edu/projects/hotnet/              |
        -------------------------------------------------------------------------------------
        |-hf/--heat_file         |Path to a tab-separated file containing a gene name in the|
        |                        |first column and the heat score for that gene in the      |
        |                        |second column of each line.                               |
        -------------------------------------------------------------------------------------
        |-pnp/                   |Path to influence matrices for permuted networks.  Include|
        |--permuted_networks_path|'NUM' in the path to be replaced with the iteration number|
        -------------------------------------------------------------------------------------

Running with only the parameters specified above will create a 'hotnet_output' directory in your
current working directory that contains 4 subdirectories each prefixed with `delta_`. Each of these
subdirectories contains results files for a different value of the delta parameter used by the
HotNet2 algorithm. The output files are:

* `components.txt`: Lists subnetworks identified as significantly altered, one per line. Genes
  within each subnetwork are separated by tabs.
* `significance.txt`: For k from 2 - 10, lists the number of subnetworks of size >= k found in the
  real data, the expected number of subnetworks of size >= k based on permuted data, and the p-value
  for the observed number of subnetworks.
* `results.json`: Contains all of the above information plus the parameters used for the run in
  JSON format to faciliate further automated processing

When using `simpleRun.py`, you may also optionally provide any or all of the parameters listed
below. If one of these parameters is not provided, it will be set to the default value shown below.

        ========================================================================================================
        | PARAMETER NAME         | DEFAULT          | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet algorithm.                                     |
        --------------------------------------------------------------------------------------------------------
        |-ccs/--min_cc_size      | 3                |Minimum size connected components that should be returned.|
        --------------------------------------------------------------------------------------------------------
        |-ms/--min_heat_score    | See description  |Minimum heat score for a gene to be eligible for inclusion|
        |                        |                  |in a returned connected component. By default, all genes  |
        |                        |                  |with positive heat scores will be included. (To include   |
        |                        |                  |genes with score zero, set min_heat_score to 0).          |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | 100              |Number of permutations that should be used for parameter  |
        |                        |                  |selection and statistical significance testing            |
        --------------------------------------------------------------------------------------------------------
        |--parallel              | Not default      |Include flag to run permutation tests in parallel. Only   |
        |                        |                  |recommended for machines with at least 8 cores.           |
        --------------------------------------------------------------------------------------------------------
        |--no-parallel           | Default          |Include flag to run permutation tests sequentially.       |
        |                        |                  |Recommended for machines with fewer than 8 cores.         |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_directory   | hotnet_output    |Output directory.                                         |
        --------------------------------------------------------------------------------------------------------

For simple runs on classic HotNet, use the `simpleRunClassic.py' Python script.  The following
parameters are required:

        =====================================================================================
        | PARAMETER NAME         | DESCRIPTION                                              |
        =====================================================================================
        |-mf/--infmat_file       |Path to .mat file containing influence matrix downloaded  |
        |                        |from http://compbio.cs.brown.edu/projects/hotnet/         |
        -------------------------------------------------------------------------------------
        |-if/--infmat_index_file |Path to gene-index mapping file downloaded from           |
        |                        |http://compbio.cs.brown.edu/projects/hotnet/              |
        -------------------------------------------------------------------------------------
        |-hf/--heat_file         |Path to a tab-separated file containing a gene name in the|
        |                        |first column and the heat score for that gene in the      |
        |                        |second column of each line.                               |
        -------------------------------------------------------------------------------------

Running with only the parameters specified above will create a 'hotnet_output' directory in your
current working directory that contains 5 subdirectories each prefixed with `delta_`. Each of these
subdirectories contains results files for a different value of the delta parameter used by the
classic HotNet algorithm.  The contents of the directories are identical to those described above
for simple runs of HotNet2 algorithm using `simpleRun.py'.

To see an example, first make sure you have downloaded the influence matrices from
[http://compbio.cs.brown.edu/projects/hotnet/](http://compbio.cs.brown.edu/projects/hotnet/)
and saved them in the `influence_matrices` directory, then run:

    python simpleRun.py @example/configs/simple.config


HotNet algorithm & code
------------------------
The HotNet algorithm, whether classic or directed, has 4 basic steps:

1. Influence matrix creation

    This step creates a matrix that defines an "influence score" for each gene pair in the network based on known gene interactions and a heat diffusion process.  Code is not provided for influence matrix creation.

2. Delta selection

    This step uses random data to select the threshold that should be used for edge weight removal in the HotNet run of step 3.
    The random data can be either:

    * Permuted influence matrices (intended for directed HotNet),  
        in which case we find the minimal delta such that |CC_{max}| <= l (by default, for l=5, 10, 15, 20),

        or

    * Permuted heat scores (intended for classic HotNet),  
        in which case we find delta that maximizes # CCs with size >= k (by default, for k=3, 4, 5)

    The file `findThreshold.py` can be used to run the delta selection step.  The required and optional parameters to the script are described in the table below.  Note that a best delta will be returned for each combination of permuted network/score and l/k.  These deltas must be manually examined to choose a single delta to use for the HotNet run in step 3.

        --------------------------------------------------------------------------------------------------------
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        --------------------------------------------------------------------------------------------------------
        |General parameters                                                                                    |
        --------------------------------------------------------------------------------------------------------
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet algorithm.                                     |
        --------------------------------------------------------------------------------------------------------
        |-p/--parallel           | None             |Include flag to run permutation tests in parallel.        |
        --------------------------------------------------------------------------------------------------------
        |-c/--classic            | None             |Include flag to run classic (instead of directed) HotNet. |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_file        | None             |Output file. If none, output will be written to stdout.   |
        --------------------------------------------------------------------------------------------------------
        |Network permutation parameters                                                                        |
        --------------------------------------------------------------------------------------------------------
        |-pnp/                   | REQUIRED         |Path to influence matrices for permuted networks. Include |
        |--permuted_networks_path|                  |##NUM## in the path path to be replaced with the iteration|
        |                        |                  |number.                                                   |
        --------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name       | Li               |Variable name of the influence matrices in the .mat files.|
        --------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file | REQUIRED         |Tab-delimited gene-index file for the influence matrices. |
        --------------------------------------------------------------------------------------------------------
        |-hf/--heat_file         | REQUIRED         |Tab-delimited heat score file.                            |
        --------------------------------------------------------------------------------------------------------
        |-k/--max_cc_sizes       | 5,10,15,20       |Max CC sizes for delta selection.                         |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted networks to use.                       |
        --------------------------------------------------------------------------------------------------------
        |Heat score permutation parameters                                                                     |
        -------------------------------------------------------------------------------------------------------
        |-mf/--infmat_file       | REQUIRED         |Path to .mat file containing influence matrix.            |
        --------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name       | Li               |Variable name of the influence matrices in the .mat files.|
        --------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file | REQUIRED         |Tab-delimited gene-index file for the influence matrix.   |
        --------------------------------------------------------------------------------------------------------
        |-hf/--heat_file         | REQUIRED         |Tab-delimited heat score file.                            |
        --------------------------------------------------------------------------------------------------------
        |-k/--max_cc_sizes       | 5,10,15,20       |Max CC sizes for delta selection.                         |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted networks to use.                       |
        --------------------------------------------------------------------------------------------------------

    _NOTE: An efficient procedure for finding the delta that maximizes the # of CCs with size >= k has not yet been implemented.
    Thus, at present, the heat permutation delta selection also finds the minimial delta such that_ |CC_{max}| <= l.

3. HotNet run

    This step performs the core HotNet algorithm: calculating a weighted graph based on influence matrix and heat score, removing edges with weight less than delta,
    and extracting the resulting connected components.

    The file `runHotnet2.py` can be used to perform the HotNet run.  The required and optional parameters to the script are described in the table below.  In addition to performing the run of the core algorithm, `runHotnet2.py` also runs the significance testing described in step 4.

        --------------------------------------------------------------------------------------------------------
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        --------------------------------------------------------------------------------------------------------
        |General parameters                                                                                    |
        --------------------------------------------------------------------------------------------------------
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet algorithm.                                     |
        --------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file | REQUIRED         |Tab-delimited gene-index file for the influence matrices. |
        --------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name       | Li               |Variable name of the influence matrices in the .mat files.|
        --------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file | REQUIRED         |Tab-delimited gene-index file for the influence matrix.   |
        --------------------------------------------------------------------------------------------------------
        |-hf/--heat_file         | REQUIRED         |Tab-delimited heat score file.                            |
        --------------------------------------------------------------------------------------------------------
        |-d/--delta              | REQUIRED         |Weight threshold for edge removal.                        |
        --------------------------------------------------------------------------------------------------------
        |-ccs/--min_cc_size      | 3                |Minimum size connected components that should be returned.|
        --------------------------------------------------------------------------------------------------------
        |Significance testing parameters                                                                       |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted networks to use. Set to 0 to skip      |
        |                        |                  |significance testing.                                     |
        --------------------------------------------------------------------------------------------------------
        |-pgf/                   | None             |Path to file containing a list of additional genes that   |
        |--permutation_genes_file|                  |can have permuted heat values assigned to them in         |
        |                        |                  |permutation tests                                         |
        --------------------------------------------------------------------------------------------------------
        |-s/--cc_start_size      | 2                |Smallest connected component size to count                |
        --------------------------------------------------------------------------------------------------------
        |-l/--cc_stop_size       | 10               |Largest connected component size to count                 |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_file        | None             |Output file. If none, output will be written to stdout.   |
        --------------------------------------------------------------------------------------------------------
        |-p/--parallel           | None             |Include flag to run significance tests in parallel.       |
        --------------------------------------------------------------------------------------------------------

4. Statistical significance testing

    This step calculates the statistical significance of obtaining the observed number of connected components of various sizes by comparing the observed number of connected components to the expected number from permuted heat score data.

    As mentioned above, statistical significance testing is included as part `runHotnet2.py` by specifying a number of permutations > 0.


Input & output
------------------------
Some or all parameters can be specified via a configuration file by passing '@<ConfigFileName>' as a command-line parameter, e.g.
`python runHotnet2.py @testConf.txt`.  Multiple configuration files can be used, and configuration files can be combined with parameters specified directly on the command line.  If the same parameter is specfied in multiple configuration files or in both a config file and on the command line, the last set value will be used.

All output is in JSON format.  It will be written to a specified output file or, if no output file is specified, to stdout.


Sample data
------------------------
Example files will be provided but are not yet available.
