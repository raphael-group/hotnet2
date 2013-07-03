HotNet2
=======================

HotNet algorithm & code
------------------------
The HotNet algorithm, whether classic or directed, has 4 basic steps:

1. Influence matrix creation

    This step creates a matrix that defines an "influence score" for each gene pair in the newtork based on known gene interactions and a heat diffusion process.  Code is not provided for influence matrix creation.

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