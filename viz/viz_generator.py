import json
import sys
sys.path.append("../")
import hotnet2 as hn
import hnio
import matplotlib.pyplot as plt

# load output from HotNet run
hn_output = json.load(open("/research/compbio/users/jeldridg/TestFiles/HotNet2_Hint_MutFreq_Output.json"))

# get lists of connected components
components = hn_output["components"]

# calculate their sizes
sizes = hn.component_sizes(components)

# verify that the calculated sizes match those in the output file (once strings in JSON keys are converted back to ints)
sizes2 = dict([(int(size), count) for size, count in hn_output["sizes"].iteritems()])
# print "Sizes from output and calculation match (expect True): %s" % (sizes == sizes2)

# get command-line parameters used for the HotNet run
parameters = hn_output["parameters"]

# print out whether this was a classic or directed HotNet run
runType = "classic" if parameters["classic"] else "directed"
# print "Was run classic or directed? (expect directed): %s" % runType 

# # load the heat scores used for the HotNet run and show them as a histogram 
heat, _ = hnio.load_heat_json(parameters["heat_file"])
# plt.hist(heat.values(), bins=100)
# plt.show()

# get parameters used for heat file generation
heat_parameters = hn_output["heat_parameters"]

# print out the number of samples with SNVs and CNAs
genes = hnio.load_genes(heat_parameters["gene_file"])
samples = hnio.load_samples(heat_parameters["sample_file"])
samples2snvs = hnio.load_snv_data(heat_parameters["snv_file"], genes, samples)
samples2cnas = hnio.load_cna_data(heat_parameters["cna_file"], genes, samples)
# print "There were %s samples with SNVs" % (len(samples2snvs))
# print "There were %s samples with CNAs" % (len(samples2cnas))

# check whether the first two genes in the largest CC interact
edges = hnio.load_ppi_edges(parameters["edge_list_file"])
index2gene = hnio.load_index(parameters["infmat_index_file"])
gene2index_mapping = {gene: index for index, gene in index2gene.items()}
gene1index = gene2index_mapping[components[0][0]]
gene2index = gene2index_mapping[components[0][1]]
genes_interact = (gene1index, gene2index) in edges or (gene2index, gene1index) in edges
# print "Do the first two genes in the largest CC interact (expect False)? %s" % (genes_interact)

# how about hte first gene and the third gene?
gene3index = gene2index_mapping[components[0][2]]
genes_interact = (gene1index, gene3index) in edges or (gene3index, gene1index) in edges
# print "Do the first and third genes in the largest CC interact (expect True)? %s" % (genes_interact)



#		END OF HOTNET CODE. IT'S VISUALIZATION TIME!

sample_list = list(samples)
gene_list = list(genes)

return_array = []

# Generating the nodes and links requires the links index into the respective nodes
# This keeps track of each node's index
index_dict = {}	

# This remembers the association of each sample to every gene that it contains
sample_dict = {sample_list[i]: {} for i in range(len(sample_list))}

# This remembers the association of each gene with every sample that contains it
gene_dict = {gene_list[i]: {} for i in range(len(gene_list))}

# Generating the association of sample: [list of genes the sample contains]
for sample in samples2snvs:
	list_of_genes = list(samples2snvs[sample])
	for i in range(len(list_of_genes)):
		cur_gene = list_of_genes[i]
		if (sample not in gene_dict[cur_gene]):
			gene_dict[cur_gene][sample] = sample
			# We populate the gene --> samples dict with the samples
		if (cur_gene not in sample_dict[sample]):
			gene_data = {	"gene": cur_gene,
							"SNV": None,
							"CNA": None,
							"locus": 100}
			sample_dict[sample][cur_gene] = gene_data
		sample_dict[sample][cur_gene]["SNV"] = True

# Still generating the association from before--this time with CNA data
for sample in samples2cnas:
	dict_of_genes = (samples2cnas[sample])
	for cur_gene in dict_of_genes:
		if (sample not in gene_dict[cur_gene]):
			gene_dict[cur_gene][sample] = sample
			# We populate the gene --> samples dict with the samples
		if (cur_gene not in sample_dict[sample]):
			gene_data = {	"gene": cur_gene,
							"SNV": None,
							"CNA": None,
							"locus": 100}
			sample_dict[sample][cur_gene] = gene_data
		sample_dict[sample][cur_gene]["CNA"] = dict_of_genes[cur_gene]

for i in range(len(components)):
	output = {	"nodes": [], 
				"links": [],
				"samples": []}
	index = 0
	cur_component = components[i]
	for j in range(len(cur_component)):
		gene1 = str(cur_component[j])
		output["nodes"].append({"gene": gene1, 
								"category": i,
								"heat": heat[gene1]})
		index_dict[gene1] = index
		index += 1

	for j in range(len(cur_component)):
		gene1 = str(cur_component[j])
		for k in range(len(cur_component)):
			gene2 = str(cur_component[k])
			if ((gene2index_mapping[gene1], gene2index_mapping[gene2]) in edges):
				output["links"].append({"source": index_dict[gene1],
										"target": index_dict[gene2],
										"present": True})
			else:
				output["links"].append({"source": index_dict[gene1],
										"target": index_dict[gene2],
										"present": False})
	for j in range(len(cur_component)):
		gene1 = str(cur_component[j])
		cur_sample_list = gene_dict[gene1].keys()
		for k in range(len(cur_sample_list)):
			cur_sample = cur_sample_list[k]
			filtered_gene_list = []
			for l in range(len(cur_component)):
				temp_cur_component = cur_component[l]
				if (temp_cur_component in sample_dict[cur_sample]):
					filtered_gene_list.append(sample_dict[cur_sample][temp_cur_component])
			output["samples"].append({	"name": cur_sample,
										"cancer": "BRCA",
										"genes": filtered_gene_list
										})

	return_array.append(output)

with open('hotnet_viz_data.json', 'w+') as outfile:
	json.dump(return_array, outfile, skipkeys=False, ensure_ascii=True, indent=1 )