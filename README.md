HotNet2
=======================

HotNet2 is an algorithm for finding significantly altered subnetworks in a large gene interaction
network. While originally developed for use with cancer mutation data, the current release also
supports any application in which scores can be assigned to genes in the network.

Requirements
------------------------

* Linux/Unix
* [Python 2.7](http://python.org/)
* [NumPy 1.6.2](http://www.numpy.org/)
* [SciPy 0.10.1](http://www.scipy.org/)
* [NetworkX 1.7](http://networkx.github.io/)

HotNet2 will likely work with additional versions of Python, NetworkX, NumPy, and SciPy, but
alternative configurations have not been tested.

Support
------------------------
For support using HotNet2, please visit the [HotNet Google Group](https://groups.google.com/forum/#!forum/hotnet-users).

Simple runs
------------------------
Follow these steps to get started running HotNet2 quickly and easily.

First, use the `makeRequiredPPRFiles.py` script to create the PPR matrix and permuted PPR matrices.
We have included configuration files for creating the required matrices for the HPRD and iRefIndex
PPI networks, which you can run as follows:

    python makeRequiredPPRFiles.py @influence_matrices/hprd.config

    python makeRequiredPPRFiles.py @influence_matrices/irefindex.config

This will create the following files in the `influence_matrices/hprd` / `influence_matrices/irefindex`
directories:

* `{hprd/iref}_index_genes`: Gene-index file for the largest component in the given network.
* `{hprd/iref}_edge_list`: Edge list file for the largest component in the given network.
* `{hprd/iref}_ppr_{alpha}.mat`: Influence matrix in MATLAB .mat format.
* `permuted`: directory containing 100 subdirectories with the above files for permuted matrices

For other networks, see the "Advanced use" section below.

Note that this step will take a long time. Fortunately, though, you only need to do it once per
interaction network you wish to use.

Once you have created the required PPR matrices, use the `simpleRun.py` Python script to run HotNet2.
You must provide the following parameters:

        =====================================================================================
        | PARAMETER NAME         | DESCRIPTION                                              |
        =====================================================================================
        |-mf/--infmat_file       |Path to .mat file containing influence matrix             |
        -------------------------------------------------------------------------------------
        |-if/--infmat_index_file |Path to gene-index mapping file                           |
        -------------------------------------------------------------------------------------
        |-hf/--heat_file         |Path to a tab-separated file containing a gene name in the|
        |                        |first column and the heat score for that gene in the      |
        |                        |second column of each line.                               |
        -------------------------------------------------------------------------------------
        |-pnp/                   |Path to influence matrices for permuted networks.  Include|
        |--permuted_networks_path|'##NUM##' in the path to be replaced with the iteration   |
        |                        |number                                                    |
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
  JSON format to faciliate further automated processing.
* `heat.json`: Input heat scores in JSON format to facilitate further automated processing.

The `simpleRun.py` script can also be used to create a web visualization of the output subnetworks.
To do so, include the `--edge_file` parameter:

        ========================================================================================================
        | PARAMETER NAME         | DEFAULT          | DESCRIPTION                                              |
        ========================================================================================================
        |-ef/--edge_file         | None             |Path to TSV file listing edges of the interaction network,|
        |                        |                  |where each row contains the indices of two genes that are |
        |                        |                  |connected in the network. This is used to create          |
        |                        |                  |subnetwork visualizations; if not provided, visualizations|
        |                        |                  |will not be made.                                         |
        --------------------------------------------------------------------------------------------------------
        |-nn/--network_name      | Network          |Display name for the interaction network.                 |
        --------------------------------------------------------------------------------------------------------

This will result in a a `viz` subdirectory of the output directory. To view the visualizations,
navigate to the `viz` directory and run `python -m SimpleHTTPServer`, then visit `http://localhost:8000`
in a browser.

When using `simpleRun.py`, you may also optionally provide any or all of the parameters listed
below. If one of these parameters is not provided, it will be set to the default value shown below.

        ========================================================================================================
        | PARAMETER NAME         | DEFAULT          | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet2 algorithm.                                    |
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

Advanced use
------------------------
The `simpleRun.py` script described above runs the entire HotNet2 pipeline in one command. For more
advanced use cases, you can also perform each step individually. In particular, you may wish to
follow the steps below when using mutation data.

The steps of the algorithm and the code provided for each step are described below.

1. ###Influence matrix creation###

    This step creates a matrix that defines an "influence score" for each gene pair in the network
    based on known gene interactions and a heat diffusion process. HotNet2 requires an influence
    matrix and a directory of permuted influence matrices as input. The release includes scripts
    for generating and permuting influence matrices.
    
    The Python script `createPPRMat.py` creates an influence matrix from a given interaction
    network. The script computes the influence matrix from the largest component in the given
    interaction network. The required and optional paramters to the script are described below.
    
        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-a/--alpha              | REQUIRED         |Page Rank dampening factor, equal to 1-beta (where beta is|
        |                        |                  |the restart probability for insulated heat diffusion      |
        |                        |                  |process).                                                 |
        --------------------------------------------------------------------------------------------------------
        |-e/--edge_list_file     | REQUIRED         |Path to a space-separated file containing a gene index in |
        |                        |                  |the first column and the gene index in the second column  |
        |                        |                  |of each line.                                             |
        --------------------------------------------------------------------------------------------------------
        |-i/--gene_index_file    | REQUIRED         |Path to a space-separated file containing a gene index in |
        |                        |                  |the first column and the corresponding gene name in the   |
        |                        |                  |second column of each line.                               |
        --------------------------------------------------------------------------------------------------------
        |-s/--start_index        | 1                |Minimum index value in the gene-index file.               |
        --------------------------------------------------------------------------------------------------------
        |--matlab                | False            |Flag whether to create the influence matrix using an      |
        |                        |                  |external call to MATLAB. Scipy is used by default. Using  |
        |                        |                  |MATLAB may lead to faster runtimes.                       |
        --------------------------------------------------------------------------------------------------------
        |-p/--output_prefix      | REQUIRED         |Prefix for each output file.                              |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_dir         | REQUIRED         |Output directory.                                         |
        --------------------------------------------------------------------------------------------------------

    The output file ares:

    * `{output_prefix}_index_genes`: Gene-index file for the largest component in the given network.
    * `{output_prefix}_edge_list`: Edge list file for the largest component in the given network.
    * `{output_prefix}_ppr_{alpha}.mat`: Influence matrix in MATLAB .mat format.

    The Python script `permuteNetwork.py` creates a directory of permuted edge lists from a given
    network. The script permutes the edges in a given interaction network -- while preserving the
    degree of each node -- by performing connected edge swaps. The required and optional paramters
    to the script are described below.
    
        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-e/--edge_list_file     | REQUIRED         |Path to a space-separated file containing a gene index in |
        |                        |                  |the first column and the gene index in the second column  |
        |                        |                  |of each line.                                             |
        --------------------------------------------------------------------------------------------------------
        |-p/--num_permutations   | 100              |Number of permuted networks to create.                    |
        --------------------------------------------------------------------------------------------------------
        |-q/--Q                  | 115              |Edge swap constant. The script attempts Q*|E| edge swaps. |
        --------------------------------------------------------------------------------------------------------
        |-s/--start_index        | 1                |Minimum index value in the gene-index file.               |
        --------------------------------------------------------------------------------------------------------
        |-p/--output_prefix      | REQUIRED         |Prefix for each output permuted edge list.                |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_dir         | REQUIRED         |Output directory.                                         |
        --------------------------------------------------------------------------------------------------------

    For convenience, we also provide the Python script `makeRequiredPPRFiles.py` for creating PPR
    matrix and an arbitrary number of permuted PPR matrices in a single step.

2. ###Heat score generation###

    This step creates a JSON file containing heat scores on each gene required in subsequent steps.
    Heat scores can either be specified directly or calculated 

    The Python script `generateHeat.py` can be used to create the JSON heat file. The required and
    optional parameters to the script are described below.

    If heat scores are specified directly, the first argument to `generateHeat.py` should be
    `scores`, e.g.:

            python generateHeat.py scores <additional_parameters>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-hf/--heat_file         | REQUIRED         |Path to a tab-separated file containing a gene name in the|
        |                        |                  |first column and the heat score for that gene in the      |
        |                        |                  |second column of each line.                               |
        --------------------------------------------------------------------------------------------------------
        |-ms/--min_heat_score    | See description  |Minimum heat score for a gene to be eligible for inclusion|
        |                        |                  |in a returned connected component. By default, all genes  |
        |                        |                  |with positive heat scores will be included. (To include   |
        |                        |                  |genes with score zero, set min_heat_score to 0).          |
        --------------------------------------------------------------------------------------------------------
        |-gff/--gene_filter_file | None             |Path to file listing genes whose heat scores should be    |
        |                        |                  |preserved. If present, heat scores for all genes not      |
        |                        |                  |listed in the file will be discarded.                     |
        --------------------------------------------------------------------------------------------------------
        |-egf/--gene_filter_file | None             |File path to which the list of genes that were excluded   |
        |--excluded_genes_output_file               |from the heat score output due to the specified filtering |
        |                        |                  |parameters should be writte, on gene per line. If no genes|
        |                        |                  |were filtered and a path is specified, the resulting file |
        |                        |                  |will be empty.                                            |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_file        | None             |Output file. If none, output will be written to stdout.   |
        --------------------------------------------------------------------------------------------------------

    If heat scores are to be calculated from mutation data, the first argument to `generateHeat.py`
    should be `mutation`, e.g.:

            python generateHeat.py mutation <additional_parameters>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |--snv_file              | REQUIRED         |Path to a tab-separated file containing single nucleotide |
        |                        |                  |variants(SNVs) where the first column of each line is a   |
        |                        |                  |sample ID and subsequent columns contain the names of     |
        |                        |                  |genes with SNVs in that sample. Lines starting with '#'   |     
        |                        |                  |will be ignored.                                          |
        --------------------------------------------------------------------------------------------------------
        |--cna_file              | REQUIRED         |Path to a tab-separated file containing copy number       |
        |                        |                  |aberrations (CNAs) where the first column of each line is |
        |                        |                  |a sample ID and subsequent columns contain gene names     |
        |                        |                  |followed by "(A)" or "(D)" indicating an amplification or |
        |                        |                  |deletion in that gene for the sample. Lines starting with |
        |                        |                  |'#' will be ignored.                                      |
        --------------------------------------------------------------------------------------------------------
        |--sample_file           | None             |Path to file listing sample IDs, one per line. Any SNVs or|
        |                        |                  |CNAs in samples not listed in this file will be ignored.  |
        |                        |                  |If HotNet2 is run with mutation permutation testing, all  |
        |                        |                  |samples in this file will be eligible for random mutations|
        |                        |                  |even if the sample did not have any mutations in the real |
        |                        |                  |data. If not provided, the set of samples is assumed to be|
        |                        |                  |all samples that are provided in the SNV or CNA data.     |
        --------------------------------------------------------------------------------------------------------
        |--gene_file             | None             |Path to file listing gene names, one per line. Mutations  |
        |                        |                  |in genes not listed in this file will be ignored. If      |
        |                        |                  |HotNet2 is run with mutation permutation testing, every   |
        |                        |                  |gene in this file will be eligible for random mutations   |
        |                        |                  |even if the gene did not have mutations in any samples in |
        |                        |                  |the original data. If not provided, the set of tested     |
        |                        |                  |genes is assumed to be all genes that have mutations in   |
        |                        |                  |either the SNV or CNA data.                               |
        --------------------------------------------------------------------------------------------------------
        |--min_freq              | 1                |The minimum number of samples in which a gene must have an|
        |                        |                  |SNV in order to be considered mutated in the heat score   |
        |                        |                  |calculation.                                              |
        --------------------------------------------------------------------------------------------------------
        |--cna_filter_threshold  | None             |Proportion of CNAs in a gene across samples that must     |
        |                        |                  |share the same CNA type in order for the CNAs to be       |
        |                        |                  |included. This must either be > .5, or the default, None, |
        |                        |                  |in which case all CNAs will be included.                  |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_file        | None             |Output file. If none, output will be written to stdout.   |
        --------------------------------------------------------------------------------------------------------

3. ###Delta selection###

    This step uses random data to select thresholds that should be used for edge weight removal in
    the HotNet2 run of step 4. The output of this step includes a lists of selected deltas for each
    permutation test.  Users must manually combine these deltas into selections for the HotNet2 run.
    Taking the median across permutation tests is recommended.
    
    The Python script `findThreshold.py` can be used to run the delta selection procedure. The
    required and optional parameters to the script are described below.
            
        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet2 algorithm.                                    |
        --------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name       | PPR              |Variable name of the influence matrices in the .mat files.|
        --------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file | REQUIRED         |Path to tab-separated file containing an index in the     |
        |                        |                  |first column and the name of the gene represented at that |
        |                        |                  |index in the second column of each line.                  |
        --------------------------------------------------------------------------------------------------------
        |-hf/--heat_file         | REQUIRED         |JSON heat score file generated via generateHeat.py        |
        --------------------------------------------------------------------------------------------------------
        |-l/--sizes              | 5, 10, 15, 20    |Values of l for selecting the smallest delta for each     |
        |                        |                  |permuted dataset such that the size of the largest CC is  |
        |                        |                  |<= l.                                                     |
        --------------------------------------------------------------------------------------------------------
        |-pnp/                   | REQUIRED         |Path to influence matrices for permuted networks.  Include|
        |--permuted_networks_path|                  |'##NUM##' in the path to be replaced with the iteration   |
        |                        |                  |number                                                    |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted data sets to generate.                 |
        --------------------------------------------------------------------------------------------------------
        |--parallel              | Not default      |Include flag to run permutation tests in parallel. Only   |
        |                        |                  |recommended for machines with at least 8 cores.           |
        --------------------------------------------------------------------------------------------------------
        |--no-parallel           | Default          |Include flag to run permutation tests sequentially.       |
        |                        |                  |Recommended for machines with fewer than 8 cores.         |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_file        | None             |Output file. If none, output will be written to stdout.   |
        --------------------------------------------------------------------------------------------------------

4. ###HotNet2 run###

    This step performs the core HotNet2 algorithm: calculating a weighted graph based on influence
    matrix and heat score, removing edges with weight less than delta, and extracting the resulting
    connected components.

    The Python script `runHotnet2.py` can be used to perform the HotNet2 run. The required and
    optional parameters to the script are described in the tables below.  In addition to performing
    the run of the core algorithm, `runHotnet2.py` also runs the significance testing described in
    step 5.

    To run HotNet2 without statistical significance testing, the final parameter to `runHotnet2.py`
    should be `none`, e.g.:

            python runHotnet2.py <additional_parameters> none

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet2 algorithm.                                    |
        --------------------------------------------------------------------------------------------------------
        |-mf/--infmat_file       | REQUIRED         |Path to .mat file containing influence matrix.            |
        --------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name       | PPR              |Variable name of the influence matrices in the .mat files.|
        --------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file | REQUIRED         |Path to tab-separated file containing an index in the     |
        |                        |                  |first column and the name of the gene represented at that |
        |                        |                  |index in the second column of each line.                  |
        --------------------------------------------------------------------------------------------------------
        |-hf/--heat_file         | REQUIRED         |JSON heat score file generated via generateHeat.py        |
        --------------------------------------------------------------------------------------------------------
        |-d/--deltas             | REQUIRED         |Weight thresholds for edge removal.                       |
        --------------------------------------------------------------------------------------------------------
        |-ccs/--min_cc_size      | 3                |Minimum size connected components that should be returned.|
        --------------------------------------------------------------------------------------------------------
        |-o/--output_directory   | REQUIRED         |Output directory. Files results.json, components.txt, and |
        |                        |                  |significance.txt will be generated                        |
        --------------------------------------------------------------------------------------------------------

5. ###Statistical significance testing###

    This step calculates the statistical significance of obtaining the observed number of connected
    components of various sizes by comparing the observed number of connected components to the
    expected number from permuted data. As with delta selection, either heat scores or mutation data
    can be be permuted to generate the random data sets.

    As mentioned above, statistical significance testing is included as part of `runHotnet2.py`. The
    additional required and optional parameters for including significance testing in the HotNet2
    run are described below.

    To run statistical significance testing with permuted heat scores, include `heat` as a parameter
    after the parameters described in step 4, then include any additional parameters for the
    permutation test, e.g.:

            python runHotNet2.py <params_from_step_4> heat <additional_significance_testing_params>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-pgf/                   | None             |Path to file containing a list of additional genes that   |
        |--permutation_genes_file|                  |can have permuted heat values assigned to them in         |
        |                        |                  |permutation tests.                                        |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted data sets to generate.                 |
        --------------------------------------------------------------------------------------------------------
        |-s/--cc_start_size      | 2                |Smallest connected component size to count                |
        --------------------------------------------------------------------------------------------------------
        |-l/--cc_stop_size       | 10               |Largest connected component size to count                 |
        --------------------------------------------------------------------------------------------------------
        |--parallel              | Not default      |Include flag to run permutation tests in parallel. Only   |
        |                        |                  |recommended for machines with at least 8 cores.           |
        --------------------------------------------------------------------------------------------------------
        |--no-parallel           | Default          |Include flag to run permutation tests sequentially.       |
        |                        |                  |Recommended for machines with fewer than 8 cores.         |
        --------------------------------------------------------------------------------------------------------

    To run statistical significance testing using pre-computed datasets, include `precomputed` as a
    parameter after the parameters described in step 4, then include any additional parameters for
    the permutation test, e.g.:

            python runHotNet2.py <params_from_step_4> mutations <additional_significance_testing_params>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-dp/--datasets_path     | REQUIRED         |Path to datasets to use for significance testing. Include |
        |                        |                  |Include ##NUM## in the path to be replaced with the       |
        |                        |                  |iteration number.                                         |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted data sets to generate.                 |
        --------------------------------------------------------------------------------------------------------
        |-s/--cc_start_size      | 2                |Smallest connected component size to count                |
        --------------------------------------------------------------------------------------------------------
        |-l/--cc_stop_size       | 10               |Largest connected component size to count                 |
        --------------------------------------------------------------------------------------------------------
        |--parallel              | Not default      |Include flag to run permutation tests in parallel. Only   |
        |                        |                  |recommended for machines with at least 8 cores.           |
        --------------------------------------------------------------------------------------------------------
        |--no-parallel           | Default          |Include flag to run permutation tests sequentially.       |
        |                        |                  |Recommended for machines with fewer than 8 cores.         |
        --------------------------------------------------------------------------------------------------------

6. Visualization

    You can visualize the subnetworks output by HotNet2 using the `makeResultsWebsite.py` script.
    It takes the following parameters:

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--results_files      | REQUIRED         |Paths to results.json files output by HotNet2. Multiple   |
        |                        |                  |file paths may be passed.                                 |
        --------------------------------------------------------------------------------------------------------
        |-ef/--edge_file         | REQUIRED         |Path to TSV file listing edges of the interaction network,|
        |                        |                  |where each row contains the indices of two genes that are |
        |                        |                  |connected in the network.                                 |
        --------------------------------------------------------------------------------------------------------
        |-nn/--network_name      | Network          |Display name for the interaction network.                 |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_directory   | REQUIRED         |Output directory.                                         |
        -------------------------------------------------------------------------------------------------------- 

    To view the resulting visualizations, navigate to the output directory and run
    `python -m SimpleHTTPServer`, then visit `http://localhost:8000` in a browser.

Passing parameters
------------------------
To avoid the need to continually pass a large number of parameters on the command line, some or all
parameters to the scripts described above can be specified via a configuration file by passing
`@<ConfigFileName>` as a command-line parameter, e.g. `python runHotnet2.py @testConf.txt`.  The
configuration file simply contains additional parameters as they would be specified on the
command line (though, unlike on the command line, line breaks are permitted between parameters in
config files). Multiple configuration files can be used, and configuration files can be combined
with parameters specified directly on the command line. If the same parameter is specified in
multiple configuration files or in both a config file and on the command line, the last set value
will be used.


Influence matrices
------------------------
In the `influence_matrices` directory, we provide edge list and gene index files for computing the
PPR matrices for the following interaction networks:

* HPRD v9 (http://www.hprd.org/)
* iRefIndex v9 (http://irefindex.org)

See above for instructions on using these files and the provided Python scripts to create the
required PPR matrices.

Citation
------------------------
If you use HotNet2 in your work, please cite:

F. Vandin, E. Upfal, and B.J. Raphael. (2011) Algorithms for Detecting Significantly Mutated
Pathways in Cancer. Journal of Computational Biology. 18(3):507-22.

F. Vandin, P. Clay, E. Upfal, and B. J. Raphael. Discovery of Mutated Subnetworks Associated with Clinical Data in Cancer. In Proc. Pacific Symposium on Biocomputing (PSB), 2012.
