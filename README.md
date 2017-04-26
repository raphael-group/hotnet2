# HotNet2 #

HotNet2 identifies subnetworks of a protein-protein interaction network with more mutations ("heat") than expected.

HotNet2 was developed by the [Raphael research group](http://compbio.cs.brown.edu/) at Brown University.

### Setup ###

#### Requirements ####

Latest tested version in parentheses:

1. Python (2.7.12)

    a. NumPy (1.12.1)

    b. SciPy (0.19.0)
    
    c. NetworkX (1.11)
    
    d. h5py (2.7.0)

2. gcc and gfortran (5.4.0)

#### Python dependencies ####
We recommend using [`virtualenv`](https://virtualenv.pypa.io/en/latest/) to install the Python requirements. After installing `virtualenv`, you can install the Python requirements for HotNet2 as follows:

    virtualenv venv
    source venv/bin/activate
    pip install -r requirements.txt

#### Compilation ###

The C and Fortran extensions are not required, but will significantly speed up HotNet2. You can compile them as follows:

    python hotnet2/setup_fortran.py build_src build_ext --inplace
    python hotnet2/setup_c.py build_src build_ext --inplace

### Usage ###

#### Recent update ####

We recently updated HotNet2 to simplify usage. This update requires some updates to your scripts and data. Consult the [previous README](https://github.com/raphael-group/hotnet2/blob/master/OLD-README.md) for comparison.

#### Data preprocessing ####

* **Heat scores**. Use `makeHeatFile.py` to create a JSON file of gene heat scores (weights). These scores can be generated from mutation data or from several other formats. `HotNet2.py` requires a path to at least one heat file as input.
* **Interaction network**. Use `makeNetworkFiles.py` to generate the network files required for running HotNet2. These include the influence matrix of your input network, and permuted networks and their corresponding influence matrices. Each of these files are in HDF5 (`.h5`) format. `HotNet2.py` requires at least one network file and path to the permuted network files as input.

See [`paper/paper_commands.sh`](https://github.com/raphael-group/hotnet2/blob/master/paper/paper_commands.sh) for examples of using the `makeHeatFile.py` and `makeNetworkFiles.py` scripts.

#### HotNet2 ####

After generating a heat file and the network files using the scripts above, use the `HotNet2.py` script to run HotNet2. The minimum arguments required for `HotNet2.py` are as follows:

    python HotNet2.py -nf <network_file> -pnp <permuted_networks_path> -hf <heat_file> -o <output_directory>

See [`paper/paper_commands.sh`](https://github.com/raphael-group/hotnet2/blob/master/paper/paper_commands.sh) for an example of using the `HotNet2.py` scripts with outputs of the `makeHeatFile.py` and `makeNetworkFiles.py` scripts.

The output of `HotNet2.py` consists of a directory containing the following:
* `{network_name}-{heat_name}/`: For each (network, heat score) pair, `HotNet2.py` outputs a directory of results. The directory contains subdirectories starting with "delta" for each delta parameter tested, each of which contain the subnetworks and statistical signifciance associated with that delta parameter.
* `consensus/`: The `consensus/` directory contains the consensus file across all networks and heat scores, 

#### Other usages ####

1. **Generate consensus from single runs**. Use `scripts/consensus.py` to generate the consensus file from the results of HotNet2 on a single network and heat score.
2. **Create dendrogram of strongly connected components**. Use `scripts/createDendrogram.py` to generate a dendrogram of strongly connected components in the HotNet2 exchanged heat graph.
3. **Permute a single network**. Use `scripts/permuteNetwork.py` to create a permuted edge list from an input network.
4. **Create influence matrix**. Use `scripts/createInfluenceMatrix.py` to create a HotNet or HotNet2 influence matrix from an input network.

### Testing ###

See [`testing/README.md`](https://github.com/raphael-group/hotnet2/blob/master/test/README.md) for testing instructions.

#### Example ####

See [`paper/paper_commands.sh`](https://github.com/raphael-group/hotnet2/blob/master/paper/paper_commands.sh) for a short but complete set of commands for reproducing the experiments in the HotNet2 paper.

### Visualization ###

We provide scripts to run an interactive web application to view the output of `HotNet2.py`, including the subnetworks in the consensus and individual runs. See [`viz/README.md`](https://github.com/raphael-group/hotnet2/blob/master/viz/README.md) and the wiki for additional instructions and details.

### Support ###

Please visit the [HotNet Google Group](https://groups.google.com/forum/#!forum/hotnet-users) to post questions and view discussions from other users about HotNet or HotNet2, or contact us through our research group's website.

### Reference ###

If you use HotNet2 in your work, please cite 

> M.D.M. Leiserson\*, F. Vandin\*, H.T. Wu, J.R. Dobson, J.V. Eldridge, J.L. Thomas, A. Papoutsaki, Y. Kim, B. Niu, M. McLellan, M.S. Lawrence, A.G. Perez, D. Tamborero, Y. Cheng, G.A. Ryslik, N. Lopez-Bigas, G. Getz, L. Ding, and B.J. Raphael. (2014) Pan-Cancer Network Analysis Identifies Combinations of Rare Somatic Mutations across Pathways and Protein Complexes. _Nature Genetics_ 47, 106–114 (2015).

If you use HotNet in your work, please cite:

> F. Vandin, E. Upfal, and B.J. Raphael. (2011) Algorithms for Detecting Significantly Mutated Pathways in Cancer. _Journal of Computational Biology_. 18(3):507-22.
> 
> F. Vandin, P. Clay, E. Upfal, and B. J. Raphael. Discovery of Mutated Subnetworks Associated with Clinical Data in Cancer. In Proc. _Pacific Symposium on Biocomputing (PSB)_, 2012.

_(* denotes equal contribution)_
