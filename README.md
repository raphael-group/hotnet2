HotNet2
=======================

HotNet2 is an algorithm for finding significantly altered subnetworks in a large gene interaction
network. While originally developed for use with cancer mutation data, the current release also
supports any application in which meaningful scores can be assigned to genes in the network.

HotNet2 vs. Classic HotNet
------------------------
This distribution contains two related algorithms: the original HotNet algorithm "classic HotNet",
and an updated version "HotNet2". For details on the two algorithms, please refer to the publications
listed at the end of this README.

In brief, HotNet2 uses a new heat diffusion kernel analogous to random walk with restart that better
captures the local topology of the interaction network as compared to the general heat diffusion
process used by classic HotNet. HotNet2 also uses an asymmetric influence score and different
permutation testing and parameter selection procedures. Although classic HotNet is included for
completeness, we recommend using HotNet2.


Requirements
------------------------

* Linux/Unix
* [Python 2.7](http://python.org/)
* [NumPy 1.6.2](http://www.numpy.org/)
* [SciPy 0.10.1](http://www.scipy.org/)
* [NetworkX 1.7](http://networkx.github.io/)
* [h5py 2.4.0](http://www.h5py.org/)
* Fortran or C compiler (optional but recommended for performance)

HotNet2 will likely work with additional versions of Python, h5py, NetworkX, NumPy, and SciPy, but
alternative configurations have not been tested.

Support
------------------------
For support using HotNet, please visit the [HotNet Google Group](https://groups.google.com/forum/#!forum/hotnet-users).

Setup
------------------------

### Compilation
For best performance, install a Fortran or C complier and run one of the following commands
(or some appropriate variation of them) prior to running HotNet for the first time:

With a Fortran compiler:

    python hotnet2/setup_fortran.py build_src build_ext --inplace

With a C compiler:

    python hotnet2/setup_c.py build_src build_ext --inplace

If you are unable to perform these steps, the code will transparently fall back to a pure Python
implementation.

### Influence matrix creation

For each gene-gene interaction network you want to use with HotNet2, you must perform a one-time
step to generate the corresponding influence matrix. Use the provided `makeRequiredPPRFiles.py`
script to create the real and permuted personalized pagerank influence matrices. We have included
configuration files for creating the required matrices for the HPRD and iRefIndex PPI networks,
which you can run as follows:

    python makeRequiredPPRFiles.py @influence_matrices/hprd.config

    python makeRequiredPPRFiles.py @influence_matrices/irefindex.config

This will create the following files in the `influence_matrices/hprd` / `influence_matrices/irefindex`
directories:

* `{hprd/iref}_index_genes`: Gene-index file for the largest component in the given network.
* `{hprd/iref}_edge_list`: Edge list file for the largest component in the given network.
* `{hprd/iref}_ppr_{beta}.h5`: Personalized page-rank influence matrix in HDF5 format.
* `permuted`: directory containing 100 subdirectories with the above files for permuted matrices

For other networks, and for creating influence matrices for use with the classic HotNet algorithm,
see the "Advanced use" section below.

Note that this step will take a long time. Fortunately, though, you only need to do it once per interaction network you wish to use. We provide 1000 permuted networks for the HINT+HI2012,
iRefIndex 9, and Multinet interaction networks for download:

* [HINT+HI2012](http://compbio-research.cs.brown.edu/software/hotnet2/permuted_networks/hint+hi2012.tar) (~500Mb)
* [iRefIndex 9](http://compbio-research.cs.brown.edu/software/hotnet2/permuted_networks/iref9.tar) (~1.2Gb)
* [Multinet](http://compbio-research.cs.brown.edu/software/hotnet2/permuted_networks/multinet.tar) (~1.4Gb)

You can use these permuted networks instead of generating your own. However, you will
still need to generate the influence matrices for each of the permuted networks. See
the "Advanced Use" section below.

Simple runs
------------------------
Once you have performed the influence matrix creation step described above, you can use the
`runHotNet2.py` script to get started running HotNet2 quickly and easily. You must provide the
following parameters:

        ========================================================================================
        | PARAMETER NAME          | DESCRIPTION                                                |
        ========================================================================================
        |-mf/--infmat_file        |Path to HDF5 (.h5) file containing influence matrix. NumPy  |
        |                         | (.np) and MATLAB files also supported.                     |
        ----------------------------------------------------------------------------------------
        |-if/--infmat_index_file  |Path to tab-separated file containing an index in the first |
        |                         |column and the name of the gene represented at that index   |
        |                         |in the second column of each line.                          |
        ----------------------------------------------------------------------------------------
        |-hf/--heat_file          |Path to heat file containing gene names and scores. This    |
        |                         |can either be a JSON file created by generateHeat.py        |
        |                         |(described below), in which case the file name must end in  |
        |                         |.json, or a  tab-separated file containing a gene name in   |
        |                         |the first column and the heat score for that gene in the    |
        |                         |second  column of each line.                                |
        ----------------------------------------------------------------------------------------
        |-pnp                     |Path to influence matrices for permuted networks. Include   |
        |--permuted_networks_path |##NUM## in the path to be replaced with the iteration       |
        |                         |number                                                      |
        ----------------------------------------------------------------------------------------

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

The `runHotNet2.py` script can also be used to create a web visualization of the output subnetworks.
To do so, include the `--edge_file` parameter, and, optionally, other visualization-related parameters:

        =============================================================================================================
        | PARAMETER NAME          | DEFAULT   | DESCRIPTION                                                         |
        =============================================================================================================
        |-ef/--edge_file          | None               |Path to TSV file listing edges of the interaction network,  |
        |                         |                    |where each row contains the indices of two genes that are   |
        |                         |                    |connected in the network. This is used to create subnetwork |
        |                         |                    |visualizations; if not provided, visualizations will not be |
        |                         |                    |made.                                                       |
        -------------------------------------------------------------------------------------------------------------
        |-dsf                     | None               |Path to a tab-separated file containing a gene name in the  |
        |--display_score_file     |                    |first column and the display score for that gene in the     |
        |                         |                    |second column of each line.                                 |
        -------------------------------------------------------------------------------------------------------------
        |-nn/--network_name       | Network            |Display name for the interaction network.                   |
        -------------------------------------------------------------------------------------------------------------

This will result in a a `viz` subdirectory of the output directory. To view the visualizations,
navigate to the `viz` directory and run `python -m SimpleHTTPServer`, then visit `http://localhost:8000`
in a browser.

When using `runHotNet2.py`, you may also optionally provide any or all of the parameters listed
below. If one of these parameters is not provided, it will be set to the default value shown below.

        =============================================================================================================
        | PARAMETER NAME          | DEFAULT            | DESCRIPTION                                                |
        =============================================================================================================
        |-r/--runname             | None               |Name of run / disease.                                      |
        -------------------------------------------------------------------------------------------------------------
        |-ccs/--min_cc_size       | 2                  |Minimum size connected components that should be returned.  |
        -------------------------------------------------------------------------------------------------------------
        |-c/--num_cores           | 1                  |Number of cores to use for running permutation tests in     |
        |                         |                    |parallel. If -1, all available cores will be used.          |
        -------------------------------------------------------------------------------------------------------------
        |-dp/--delta_permutations | 100                |Number of permutations to be used for delta parameter       |
        |                         |                    |selection.                                                  |
        -------------------------------------------------------------------------------------------------------------
        |-sp                      | 100                |Number of permutations to be used for statistical           |
        |--significance_permutatio|                    |significance testing.                                       |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_directory    | hotnet_output      |Output directory. Files results.json, components.txt, and   |
        |                         |                    |significance.txt will be generated in subdirectories for    |
        |                         |                    |each delta.                                                 |
        -------------------------------------------------------------------------------------------------------------

For simple runs on classic HotNet, use the `runClassicHotNet.py` Python script.  The following
parameters are required:

        ========================================================================================
        | PARAMETER NAME          | DESCRIPTION                                                |
        ========================================================================================
        |-mf/--infmat_file        |Path to .mat file containing influence matrix               |
        ----------------------------------------------------------------------------------------
        |-if/--infmat_index_file  |Path to tab-separated file containing an index in the first |
        |                         |column and the name of the gene represented at that index   |
        |                         |in the second column of each line.                          |
        ----------------------------------------------------------------------------------------
        |-hf/--heat_file          |Path to heat file containing gene names and scores. This    |
        |                         |can either be a JSON file created by generateHeat.py        |
        |                         |(described below), in which case the file name must end in  |
        |                         |.json, or a  tab-separated file containing a gene name in   |
        |                         |the first column and the heat score for that gene in the    |
        |                         |second  column of each line.                                |
        ----------------------------------------------------------------------------------------
        |-pnp                     |Path to influence matrices for permuted networks. Include   |
        |--permuted_networks_path |##NUM## in the path to be replaced with the iteration       |
        |                         |number                                                      |
        ----------------------------------------------------------------------------------------

Running with only the parameters specified above will create a 'hotnet_output' directory in your
current working directory that contains 5 subdirectories each prefixed with `delta_`. Each of these
subdirectories contains results files for a different value of the delta parameter used by the
classic HotNet algorithm.  The contents of the directories are identical to those described above
for simple runs of HotNet2 algorithm using `runHotNet2.py`. Similarily, the `runClassicHotNet.py`
script accepts the same optional parameters as the `runHotNet2.py` script.


Advanced use
------------------------
The `runHotNet2.py` and `runClassicHotNet.py` scripts described above run the entire HotNet
pipeline in one command. For more advanced use cases, you can also perform each step individually.
In particular, you may wish to follow the steps below when using mutation data.  Note that the
code described below is used for both HotNet2 and classic HotNet; in all cases, passing a
`--classic` command line flag produces the classic HotNet behavior.

The steps of the algorithm and the code provided for each step are described below.

1. ###Influence matrix creation###

    This step creates a matrix that defines an "influence score" for each gene pair in the network
    based on known gene interactions and a heat diffusion process. As mentioned above, this step
    only needs to be performed *once* for each gene interaction network.

    ####HotNet2 influence matrices####
    HotNet2 requires an influence matrix and a directory of permuted influence matrices as input.
    The script `makeRequiredPPRFiles.py` can be used to create the real and permuted matrices.
    The required and optional parameters to the script are as follows:

        =============================================================================================================
        | PARAMETER NAME          | REQUIRED/DEFAULT   | DESCRIPTION                                                |
        =============================================================================================================
        |-e/--edgelist_file       | REQUIRED           |Path to TSV file listing edges of the interaction network,  |
        |                         |                    |where each row contains the indices of two genes that are   |
        |                         |                    |connected in the network.                                   |
        -------------------------------------------------------------------------------------------------------------
        |-i/--gene_index_file     | REQUIRED           |Path to tab-separated file containing an index in the first |
        |                         |                    |column and the name of the gene represented at that index   |
        |                         |                    |in the second column of each line.                          |
        -------------------------------------------------------------------------------------------------------------
        |-p/--prefix              | REQUIRED           |Output prefix.                                              |
        -------------------------------------------------------------------------------------------------------------
        |-is                      | 1                  |Minimum index in the index file.                            |
        |--index_file_start_index |                    |                                                            |
        -------------------------------------------------------------------------------------------------------------
        |-b/--beta                | REQUIRED           |Beta is the restart probability for the insulated heat      |
        |                         |                    |diffusion process.                                          |
        -------------------------------------------------------------------------------------------------------------
        |-q/--Q                   | 115                |Edge swap constant. The script will attempt Q*|E| edge      |
        |                         |                    |swaps                                                       |
        -------------------------------------------------------------------------------------------------------------
        |-ps                      | 1                  |Index at which to start permutation file names.             |
        |--permutation_start_index|                    |                                                            |
        -------------------------------------------------------------------------------------------------------------
        |-n/--num_permutations    | 100                |Number of permuted networks to create.                      |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_dir          | REQUIRED           |Output directory.                                           |
        -------------------------------------------------------------------------------------------------------------

    If desired, the scripts `bin/createPPRMat.py` and `bin/permuteNetwork.py` can be used to
    perform the individual steps of creating influence matrices and permuting edge lists,
    respectively.

    ####Classic HotNet influence matrices####
    Classic HotNet does not use permuted influence matrices, and thus requires only the single real
    influence matrix. This can be created using the `bin/createClassicInfmat.py` script. The required
    and optional parameters to the script are as follows:

        =============================================================================================================
        | PARAMETER NAME          | REQUIRED/DEFAULT   | DESCRIPTION                                                |
        =============================================================================================================
        |-e/--edgelist_file       | REQUIRED           |Path to TSV file listing edges of the interaction network,  |
        |                         |                    |where each row contains the indices of two genes that are   |
        |                         |                    |connected in the network.                                   |
        -------------------------------------------------------------------------------------------------------------
        |-i/--gene_index_file     | REQUIRED           |Path to tab-separated file containing an index in the first |
        |                         |                    |column and the name of the gene represented at that index   |
        |                         |                    |in the second column of each line.                          |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_dir          | REQUIRED           |Path to output directory.                                   |
        -------------------------------------------------------------------------------------------------------------
        |-p/--prefix              | REQUIRED           |Output prefix.                                              |
        -------------------------------------------------------------------------------------------------------------
        |-s/--start_index         | 1                  |Minimum index in the index file.                            |
        -------------------------------------------------------------------------------------------------------------
        |-t/--time                | REQUIRED           |Diffusion time.                                             |
        -------------------------------------------------------------------------------------------------------------


2. ###Heat score generation###

    This step creates a JSON file containing heat scores on each gene required in subsequent steps.
    Heat scores can either be specified directly or calculated from mutation, MutSig, MuSic, or
    Oncodrive data.

    The Python script `generateHeat.py` can be used to create the JSON heat file. The required and
    optional parameters to the script are described below.

    If heat scores are specified directly, the first argument to `generateHeat.py` should be
    `scores`, e.g.:

            python generateHeat.py scores <additional_parameters>

        =============================================================================================================
        | PARAMETER NAME          | REQUIRED/DEFAULT   | DESCRIPTION                                                |
        =============================================================================================================
        |-hf/--heat_file          | REQUIRED           |Path to a tab-separated file containing a gene name in the  |
        |                         |                    |first column and the heat score for that gene in the second |
        |                         |                    |column of each line.                                        |
        -------------------------------------------------------------------------------------------------------------
        |-ms/--min_heat_score     | 0                  |Minimum heat score for genes to have their original heat    |
        |                         |                    |score in the resulting output file. Genes with score below  |
        |                         |                    |this value will be assigned score 0.                        |
        -------------------------------------------------------------------------------------------------------------
        |-gff/--gene_filter_file  | None               |Path to file listing genes whose heat scores should be      |
        |                         |                    |preserved, one per line. If present, all other heat scores  |
        |                         |                    |will be discarded.                                          |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_file         | None               |Output file. If none given, output will be written to       |
        |                         |                    |stdout.                                                     |
        -------------------------------------------------------------------------------------------------------------

    If heat scores are to be calculated from mutation data, the first argument to `generateHeat.py`
    should be `mutation`, e.g.:

            python generateHeat.py mutation <additional_parameters>

        =============================================================================================================
        | PARAMETER NAME          | REQUIRED/DEFAULT   | DESCRIPTION                                                |
        =============================================================================================================
        |--snv_file               | REQUIRED           |Path to a tab-separated file containing SNVs where the      |
        |                         |                    |first column of each line is a sample ID and subsequent     |
        |                         |                    |columns contain the names of genes with SNVs in that        |
        |                         |                    |sample. Lines starting with "#" will be ignored.            |
        -------------------------------------------------------------------------------------------------------------
        |--cna_file               | None               |Path to a tab-separated file containing CNAs where the      |
        |                         |                    |first column of each line is a sample ID and subsequent     |
        |                         |                    |columns contain gene names followed by "(A)" or "(D)"       |
        |                         |                    |indicating an amplification or deletion in that gene for    |
        |                         |                    |the sample. Lines starting with "#" will be ignored.        |
        -------------------------------------------------------------------------------------------------------------
        |--sample_file            | None               |File listing samples. Any SNVs or CNAs in samples not       |
        |                         |                    |listed in this file will be ignored. If HotNet is run with  |
        |                         |                    |mutation permutation testing, all samples in this file will |
        |                         |                    |be eligible for random mutations even if the sample did not |
        |                         |                    |have any mutations in the real data. If not provided, the   |
        |                         |                    |set of samples is assumed to be all samples that are        |
        |                         |                    |provided in the SNV or CNA data.                            |
        -------------------------------------------------------------------------------------------------------------
        |--sample_type_file       | None               |File listing type (e.g. cancer, datasets, etc.) of samples  |
        |                         |                    |(see --sample_file). Each line is a space-separated row     |
        |                         |                    |listing one sample and its type. The sample types are used  |
        |                         |                    |for creating the HotNet(2) web output.                      |
        -------------------------------------------------------------------------------------------------------------
        |--gene_file              | None               |File listing tested genes. SNVs or CNAs in genes not listed |
        |                         |                    |in this file will be ignored. If HotNet is run with         |
        |                         |                    |mutation permutation testing, every gene in this file will  |
        |                         |                    |be eligible for random mutations even if the gene did not   |
        |                         |                    |have mutations in any samples in the original data. If not  |
        |                         |                    |provided, the set of tested genes is assumed to be all      |
        |                         |                    |genes that have mutations in either the SNV or CNA data.    |
        -------------------------------------------------------------------------------------------------------------
        |--min_freq               | 1                  |The minimum number of samples in which a gene must have an  |
        |                         |                    |SNV to be considered mutated in the heat score calculation. |
        -------------------------------------------------------------------------------------------------------------
        |--cna_filter_threshold   | None               |Proportion of CNAs in a gene across samples that must share |
        |                         |                    |the same CNA type in order for the CNAs to be included.     |
        |                         |                    |This must either be > .5, or the default, None, in which    |
        |                         |                    |case all CNAs will be included.                             |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_file         | None               |Output file. If none given, output will be written to       |
        |                         |                    |stdout.                                                     |
        -------------------------------------------------------------------------------------------------------------

    If heat scores are to be calculated from MutSig data, the first argument to `generateHeat.py`
    should be `mutsig`, e.g.:

            python generateHeat.py mutsig <additional_parameters>

        =============================================================================================================
        | PARAMETER NAME          | REQUIRED/DEFAULT   | DESCRIPTION                                                |
        =============================================================================================================
        |--mutsig_score_file      | REQUIRED           |MutSig score file (gene to q-value).                        |
        -------------------------------------------------------------------------------------------------------------
        |--threshold              | 1.0                |Maximum q-value threshold.                                  |
        -------------------------------------------------------------------------------------------------------------
        |--gene_filter_file       | None               |File listing genes whose heat scores should be preserved.   |
        |                         |                    |If present, all other heat scores will be discarded.        |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_file         | None               |Output file. If none given, output will be written to       |
        |                         |                    |stdout.                                                     |
        -------------------------------------------------------------------------------------------------------------

    If heat scores are to be calculated from MuSiC data, the first argument to `generateHeat.py`
    should be `music`, e.g.:

            python generateHeat.py music <additional_parameters>

        =============================================================================================================
        | PARAMETER NAME          | REQUIRED/DEFAULT   | DESCRIPTION                                                |
        =============================================================================================================
        |--music_score_file       | REQUIRED           |MuSiC score file (gene to q-value).                         |
        -------------------------------------------------------------------------------------------------------------
        |--threshold              | 1.0                |Maximum q-value threshold.                                  |
        -------------------------------------------------------------------------------------------------------------
        |--max_heat               | 15                 |Max heat                                                    |
        -------------------------------------------------------------------------------------------------------------
        |--gene_filter_file       | None               |File listing genes whose heat scores should be preserved.   |
        |                         |                    |If present, all other heat scores will be discarded.        |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_file         | None               |Output file. If none given, output will be written to       |
        |                         |                    |stdout.                                                     |
        -------------------------------------------------------------------------------------------------------------

    If heat scores are to be calculated from Oncodrive data, the first argument to `generateHeat.py`
    should be `oncodrive`, e.g.:

            python generateHeat.py oncodrive <additional_parameters>

        =============================================================================================================
        | PARAMETER NAME          | REQUIRED/DEFAULT   | DESCRIPTION                                                |
        =============================================================================================================
        |--fm_scores              | REQUIRED           |Oncodrive-FM scores (gene to q-value).                      |
        -------------------------------------------------------------------------------------------------------------
        |--cis_amp_scores         | REQUIRED           |Oncodrive-CIS scores (gene to q-value); amplifications      |
        |                         |                    |only.                                                       |
        -------------------------------------------------------------------------------------------------------------
        |--cis                    | False              |Flag whether to include Oncodrive-CIS scores when           |
        |                         |                    |generating the Oncodrive heat file.                         |
        -------------------------------------------------------------------------------------------------------------
        |--cis_del_scores         | REQUIRED           |Oncodrive-CIS scores (gene to q-value); deletions only.     |
        -------------------------------------------------------------------------------------------------------------
        |--fm_threshold           | 0.2                |Maximum Oncodrive-FM q-value threshold                      |
        -------------------------------------------------------------------------------------------------------------
        |--cis_threshold          | 0.2                |Maximum Oncodrive-CIS q-value threshold                     |
        -------------------------------------------------------------------------------------------------------------
        |--gene_filter_file       | None               |File listing genes whose heat scores should be preserved.   |
        |                         |                    |If present, all other heat scores will be discarded.        |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_file         | None               |Output file. If none given, output will be written to       |
        |                         |                    |stdout.                                                     |
        -------------------------------------------------------------------------------------------------------------


3. ###Delta selection###

    This step uses random data to select thresholds that should be used for edge weight removal in
    the HotNet run of step 4. The output of this step includes a lists of selected deltas for each
    permutation test.  Users must manually combine these deltas into selections for the HotNet run.
    Taking the median across permutation tests is recommended.

    The random data can be either:

    * Permuted networks (recommended for HotNet2),

    * Permuted heat scores (recommended for classic HotNet with arbitrary score types),

      or

    * Permuted mutation data (recommended for classic HotNet when mutation data is used to generate scores)

    The Python script `bin/findThreshold.py` can be used to run the delta selection procedure. The
    required and optional parameters to the script are described below.

    For the permuted networks test, precomputed network permutations are required as input. In this
    case, the first parameter to `findThreshold.py` should be `network`, e.g.:

            python bin/findThreshold.py network <additional_parameters>

        =============================================================================================================
        | PARAMETER NAME          | REQUIRED/DEFAULT   | DESCRIPTION                                                |
        =============================================================================================================
        |-r/--runname             | None               |Name of run / disease.                                      |
        -------------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name        | PPR                |Variable name of the influence matrices in the .mat files   |
        -------------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file  | REQUIRED           |Path to tab-separated file containing an index in the first |
        |                         |                    |column and the name of the gene represented at that index   |
        |                         |                    |in the second column of each line.                          |
        -------------------------------------------------------------------------------------------------------------
        |-hf/--heat_file          | REQUIRED           |JSON heat score file generated via generateHeat.py          |
        -------------------------------------------------------------------------------------------------------------
        |-n/--num_permutations    | REQUIRED           |Number of permuted data sets to generate                    |
        -------------------------------------------------------------------------------------------------------------
        |-s/--test_statistic      | max_cc_size        |If max_cc_size, select smallest delta for each permuted     |
        |                         |                    |dataset such that the size of the largest CC is <= l. If    |
        |                         |                    |num_ccsselect for each permuted dataset the delta that      |
        |                         |                    |maximizes the number of CCs of size >= l.                   |
        -------------------------------------------------------------------------------------------------------------
        |-l/--sizes               | 5, 10, 15, 20      |See test_statistic. For test_statistic 'num_ccs', default   |
        |                         |                    |is 3.                                                       |
        -------------------------------------------------------------------------------------------------------------
        |-c/--num_cores           | 1                  |Number of cores to use for running permutation tests in     |
        |                         |                    |parallel. If -1, all available cores will be used.          |
        -------------------------------------------------------------------------------------------------------------
        |--classic                | False              |Run classic (instead of directed) HotNet.                   |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_file         | None               |Output file. If none given, output will be written to       |
        |                         |                    |stdout.                                                     |
        -------------------------------------------------------------------------------------------------------------
        |-pnp                     | REQUIRED           |Path to influence matrices for permuted networks. Include   |
        |--permuted_networks_path |                    |##NUM## in the path to be replaced with the iteration       |
        |                         |                    |number                                                      |
        -------------------------------------------------------------------------------------------------------------


    For the permuted heat scores test, heat scores will simply be shuffled among genes. In this
    case, the first parameter to `findTheshold.py` should be `heat`, e.g.:

            python bin/findThreshold.py heat <additional_parameters>

        =============================================================================================================
        | PARAMETER NAME          | REQUIRED/DEFAULT   | DESCRIPTION                                                |
        =============================================================================================================
        |-r/--runname             | None               |Name of run / disease.                                      |
        -------------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name        | PPR                |Variable name of the influence matrices in the .mat files   |
        -------------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file  | REQUIRED           |Path to tab-separated file containing an index in the first |
        |                         |                    |column and the name of the gene represented at that index   |
        |                         |                    |in the second column of each line.                          |
        -------------------------------------------------------------------------------------------------------------
        |-hf/--heat_file          | REQUIRED           |JSON heat score file generated via generateHeat.py          |
        -------------------------------------------------------------------------------------------------------------
        |-n/--num_permutations    | REQUIRED           |Number of permuted data sets to generate                    |
        -------------------------------------------------------------------------------------------------------------
        |-s/--test_statistic      | max_cc_size        |If max_cc_size, select smallest delta for each permuted     |
        |                         |                    |dataset such that the size of the largest CC is <= l. If    |
        |                         |                    |num_ccsselect for each permuted dataset the delta that      |
        |                         |                    |maximizes the number of CCs of size >= l.                   |
        -------------------------------------------------------------------------------------------------------------
        |-l/--sizes               | 5, 10, 15, 20      |See test_statistic. For test_statistic 'num_ccs', default   |
        |                         |                    |is 3.                                                       |
        -------------------------------------------------------------------------------------------------------------
        |-c/--num_cores           | 1                  |Number of cores to use for running permutation tests in     |
        |                         |                    |parallel. If -1, all available cores will be used.          |
        -------------------------------------------------------------------------------------------------------------
        |--classic                | False              |Run classic (instead of directed) HotNet.                   |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_file         | None               |Output file. If none given, output will be written to       |
        |                         |                    |stdout.                                                     |
        -------------------------------------------------------------------------------------------------------------
        |-mf/--infmat_file        | REQUIRED           |Path to .mat file containing influence matrix               |
        -------------------------------------------------------------------------------------------------------------
        |-pgf                     | None               |Path to file containing a list of additional genes that can |
        |--permutation_genes_file |                    |have permuted heat values assigned to them in permutation   |
        |                         |                    |tests                                                       |
        -------------------------------------------------------------------------------------------------------------


    For the permuted mutation data test, SNVs will be randomly generated in genes according to a
    specified background mutation rate, and each CNA block will be randomly placed on an equally-
    sized group of genes in the same chromosome. In this case, the first parameter to `findTheshold.py`
    should be `mutations`, e.g.:

        	python bin/findThreshold.py mutations <additional_parameters>

        =============================================================================================================
        | PARAMETER NAME          | REQUIRED/DEFAULT   | DESCRIPTION                                                |
        =============================================================================================================
        |-r/--runname             | None               |Name of run / disease.                                      |
        -------------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name        | PPR                |Variable name of the influence matrices in the .mat files   |
        -------------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file  | REQUIRED           |Path to tab-separated file containing an index in the first |
        |                         |                    |column and the name of the gene represented at that index   |
        |                         |                    |in the second column of each line.                          |
        -------------------------------------------------------------------------------------------------------------
        |-hf/--heat_file          | REQUIRED           |JSON heat score file generated via generateHeat.py          |
        -------------------------------------------------------------------------------------------------------------
        |-n/--num_permutations    | REQUIRED           |Number of permuted data sets to generate                    |
        -------------------------------------------------------------------------------------------------------------
        |-s/--test_statistic      | max_cc_size        |If max_cc_size, select smallest delta for each permuted     |
        |                         |                    |dataset such that the size of the largest CC is <= l. If    |
        |                         |                    |num_ccsselect for each permuted dataset the delta that      |
        |                         |                    |maximizes the number of CCs of size >= l.                   |
        -------------------------------------------------------------------------------------------------------------
        |-l/--sizes               | 5, 10, 15, 20      |See test_statistic. For test_statistic 'num_ccs', default   |
        |                         |                    |is 3.                                                       |
        -------------------------------------------------------------------------------------------------------------
        |-c/--num_cores           | 1                  |Number of cores to use for running permutation tests in     |
        |                         |                    |parallel. If -1, all available cores will be used.          |
        -------------------------------------------------------------------------------------------------------------
        |--classic                | False              |Run classic (instead of directed) HotNet.                   |
        -------------------------------------------------------------------------------------------------------------
        |-o/--output_file         | None               |Output file. If none given, output will be written to       |
        |                         |                    |stdout.                                                     |
        -------------------------------------------------------------------------------------------------------------
        |-mf/--infmat_file        | REQUIRED           |Path to .mat file containing influence matrix               |
        -------------------------------------------------------------------------------------------------------------
        |-glf/--gene_length_file  | REQUIRED           |Path to tab-separated file containing gene names in the     |
        |                         |                    |first column and the length of the gene in base pairs in    |
        |                         |                    |the second column                                           |
        -------------------------------------------------------------------------------------------------------------
        |-gof/--gene_order_file   | REQUIRED           |Path to file containing tab-separated lists of genes on     |
        |                         |                    |each chromosme, in order of their position on the           |
        |                         |                    |chromosome, one chromosome per line                         |
        -------------------------------------------------------------------------------------------------------------
        |-b/--bmr                 | REQUIRED           |Default background mutation rate                            |
        -------------------------------------------------------------------------------------------------------------
        |-bf/--bmr_file           | None               |File listing gene-specific BMRs. If none, the default BMR   |
        |                         |                    |will be used for all genes.                                 |
        -------------------------------------------------------------------------------------------------------------


4. ###HotNet run###

    This step performs the core HotNet algorithm: calculating a weighted graph based on influence
    matrix and heat score, removing edges with weight less than delta, and extracting the resulting
    connected components.

    The Python script `bin/findComponents.py` can be used to perform the HotNet run.  The required and
    optional parameters to the script are described in the tables below.  In addition to performing
    the run of the core algorithm, `findComponents.py` also runs the significance testing described in
    step 5.

    To run HotNet without statistical significance testing, the final parameter to `findComponents.py`
    should be `none`, e.g.:

            python bin/findComponents.py <additional_parameters> none

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-r/--runname            | None             |Name of run/disease. This is only for the user's          |
        |                        |                  |record-keeping convenience; the parameter is not used by  |
        |                        |                  |the HotNet algorithm.                                     |
        --------------------------------------------------------------------------------------------------------
        |-mf/--infmat_file       | REQUIRED         |Path to .mat file containing influence matrix.            |
        --------------------------------------------------------------------------------------------------------
        |-mn/--infmat_name       | PPR              |Variable name of the influence matrix in the .mat file.   |
        --------------------------------------------------------------------------------------------------------
        |-if/--infmat_index_file | REQUIRED         |Path to tab-separated file containing an index in the     |
        |                        |                  |first column and the name of the gene represented at that |
        |                        |                  |index in the second column of each line.                  |
        --------------------------------------------------------------------------------------------------------
        |-hf/--heat_file         | REQUIRED         |JSON heat score file generated via generateHeat.py        |
        --------------------------------------------------------------------------------------------------------
        |-d/--deltas             | REQUIRED         |Weight thresholds for edge removal.                       |
        --------------------------------------------------------------------------------------------------------
        |-ccs/--min_cc_size      | 2                |Minimum size connected components that should be returned.|
        --------------------------------------------------------------------------------------------------------
        |-o/--output_directory   | REQUIRED         |Output directory. Files results.json, components.txt, and |
        |                        |                  |significance.txt will be generated in subdirectories for  |
        |                        |                  |each delta.                                               |
        --------------------------------------------------------------------------------------------------------
        |-c/--classic            | None             |Run classic HotNet (rather than HotNet2).                 |
        --------------------------------------------------------------------------------------------------------

5. ###Statistical significance testing###

    This step calculates the statistical significance of obtaining the observed number of connected
    components of various sizes by comparing the observed number of connected components to the
    expected number from permuted data. As with delta selection, either heat scores or mutation data
    can be be permuted to generate the random data sets.

    As mentioned above, statistical significance testing is included as part of `findComponents.py`. The
    additional required and optional parameters for including significance testing in the HotNet
    run are described below.

    To run statistical significance testing with permuted heat scores, include `heat` as a parameter
    after the parameters described in step 4, then include any additional parameters for the
    permutation test, e.g.:

            python bin/findComponents.py <params_from_step_4> heat <additional_significance_testing_params>

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
        |-c/--num_cores          | 1                |Number of cores to use for running permutation tests in   |
        |                        |                  |parallel. If -1, all available cores will be used.        |
        --------------------------------------------------------------------------------------------------------

    To run statistical significance testing with permuted mutation data, include `mutations` as a
    parameter after the parameters described in step 4, then include any additional parameters for
    the permutation test, e.g.:

            python runHotNet.py <params_from_step_4> mutations <additional_significance_testing_params>

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-glf/--gene_length_file | REQUIRED         |Path to tab-separated file containing gene names in the   |
        |                        |                  |first column and the length of the gene in base pairs in  |
        |                        |                  |the second column                                         |
        --------------------------------------------------------------------------------------------------------
        |-gof/--gene_order_file  | REQUIRED         |Path to file containing tab-separated lists of genes on   |
        |                        |                  |each chromosme, in order of their position on the         |
        |                        |                  |chromosome, one chromosome per line                       |
        --------------------------------------------------------------------------------------------------------
        |-b/--bmr                | REQUIRED         |Default background mutation rate                          |
        --------------------------------------------------------------------------------------------------------
        |-bf/--bmr_file          | None             |File listing gene-specific BMRs. If none, the default BMR |
        |                        |                  |will be used for all genes.                               |
        --------------------------------------------------------------------------------------------------------
        |-n/--num_permutations   | REQUIRED         |Number of permuted data sets to generate.                 |
        --------------------------------------------------------------------------------------------------------
        |-s/--cc_start_size      | 2                |Smallest connected component size to count                |
        --------------------------------------------------------------------------------------------------------
        |-l/--cc_stop_size       | 10               |Largest connected component size to count                 |
        --------------------------------------------------------------------------------------------------------
        |-c/--num_cores          | 1                |Number of cores to use for running permutation tests in   |
        |                        |                  |parallel. If -1, all available cores will be used.        |
        --------------------------------------------------------------------------------------------------------

6. ###Visualization###

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
        |-dsf/                   | None             |Path to a tab-separated file containing a gene name in the|
        |--display_score_file    |                  |first column and the display score for that gene in the   |
        |                        |                  |second column of each line.                               |
        --------------------------------------------------------------------------------------------------------
        |-nn/--network_name      | Network          |Display name for the interaction network.                 |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_directory   | REQUIRED         |Output directory.                                         |
        --------------------------------------------------------------------------------------------------------

    To view the resulting visualizations, navigate to the output directory and run
    `python -m SimpleHTTPServer`, then visit `http://localhost:8000` in a browser.

7. ###Consensus procedure###

    You can produce consensus output from multiple HotNet2 runs (e.g. on different protein-protein
    interaction networks) using the procedure described in the HotNet2 paper (see References below). To do
    so, use the `identifyConsensus.py` script, which takes the following parameters:

        ========================================================================================================
        | PARAMETER NAME         | REQUIRED/DEFAULT | DESCRIPTION                                              |
        ========================================================================================================
        |-d/--directories        | None             |Paths to directories with results.json files output by    |
        |                        |                  |HotNet2. The consensus will choose one results.json file  |
        |                        |                  |per directory. Only one of directories and files may be   |
        |                        |                  |set.                                                      |
        --------------------------------------------------------------------------------------------------------
        |-f/--files              | None             |Paths to (multiple) results.json files output by HotNet2. |
        |                        |                  |Only one of directories and files may be set.             |
        --------------------------------------------------------------------------------------------------------
        |-n/--networks           | REQUIRED         |List of network names, one per directory or result file.  |
        --------------------------------------------------------------------------------------------------------
        |-p/--p_value_threshold  | 0.01             |Threshold for p-values when selecting results files from  |
        |                        |                  |directories.                                              |
        --------------------------------------------------------------------------------------------------------
        |-m/--min_cc_size        | 2                |Restrict consensus to subnetworks of at least this size.  |
        --------------------------------------------------------------------------------------------------------
        |-o/--output_file        | REQUIRED         |Output file. A json extension creates a JSON file.        |
        --------------------------------------------------------------------------------------------------------

    The output file is a table where each line is a consensus subnetwork of the form `"[Core]\tExtended"`.
    The core nodes are those found in every network, and extended nodes are those found in at least one network.
    The groups of nodes are `", "`-delimited. Note that you can provide multiple results files per network
    (i.e. using different heat scores), so finding core nodes are not necessarily found in every results file.
    The output file also includes a list of "linker" genes, which are genes connected to multiple consensus subnetworks.


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


Example data
------------------------
In the `example` directory, we provide example heat, interaction network, and mutation data files,
as well as example HotNet2 configuration files. Please consult `example/README.md` for further
details.


Citation
------------------------
If you use HotNet2 in your work, please cite (\* denotes equal contribution):

M.D.M. Leiserson\*, F. Vandin\*, H.T. Wu, J.R. Dobson, J.V. Eldridge, J.L. Thomas, A. Papoutsaki,
Y. Kim, B. Niu, M. McLellan, M.S. Lawrence, A.G. Perez, D. Tamborero, Y. Cheng, G.A. Ryslik,
N. Lopez-Bigas, G. Getz, L. Ding, and B.J. Raphael.  (2014) Pan-Cancer Network Analysis Identifies
Combinations of Rare Somatic Mutations across Pathways and Protein Complexes. Nature Genetics **47**, 106–114 (2015).

If you use HotNet in your work, please cite:

F. Vandin, E. Upfal, and B.J. Raphael. (2011) Algorithms for Detecting Significantly Mutated
Pathways in Cancer. Journal of Computational Biology. 18(3):507-22.

F. Vandin, P. Clay, E. Upfal, and B. J. Raphael. Discovery of Mutated Subnetworks Associated with Clinical Data in Cancer. In Proc. Pacific Symposium on Biocomputing (PSB), 2012.
