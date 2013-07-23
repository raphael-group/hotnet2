import json
import sys
sys.path.append("../")
import hotnet2 as hn
import hnio
import matplotlib.pyplot as plt

# load output from HotNet run
hn_output = json.load(open("/research/compbio/users/jeldridg/PanCancer/Experiments/2013-07-22/FromMyCode/hint-freq.json"))

# get lists of connected components
components = hn_output["components"]

# calculate their sizes
sizes = hn.component_sizes(components)

# verify that the calculated sizes match those in the output file (once strings in JSON keys are converted back to ints)
sizes2 = dict([(int(size), count) for size, count in hn_output["sizes"].iteritems()])
print "Sizes from output and calculation match (expect True): %s" % (sizes == sizes2)

# get command-line parameters used for the HotNet run
parameters = hn_output["parameters"]

# print out whether this was a classic or directed HotNet run
runType = "classic" if parameters["classic"] else "directed"
print "Was run classic or directed? (expect directed): %s" % runType 

# load the heat scores used for the HotNet run and show them as a histogram 
heat, _ = hnio.load_heat_json(parameters["heat_file"])
plt.hist(heat.values(), bins=100)
plt.show()

# get parameters used for heat file generation
heat_parameters = hn_output["heat_parameters"]

# print out the number of samples with SNVs and CNAs
genes = hnio.load_genes(heat_parameters["gene_file"])
samples = hnio.load_samples(heat_parameters["sample_file"])
samples2snvs = hnio.load_snv_data(heat_parameters["snv_file"], genes, samples)
samples2cnas = hnio.load_cna_data(heat_parameters["cna_file"], genes, samples)
print "There were %s samples with SNVs" % (len(samples2snvs))
print "There were %s samples with CNAs" % (len(samples2cnas))
