import json
import sys
sys.path.append("../")
import hotnet2 as hn
import hnio
import matplotlib.pyplot as plt
import random

# load output from HotNet run
# hn_output = json.load(open("run2.json"))
hn_output = json.load(open("/research/compbio/users/jeldridg/TestFiles/HotNet2_Hint_MutFreq_Output.json"))

# get lists of connected components
components = hn_output["components"]

# calculate their sizes
# sizes = hn.component_sizes(components)

# verify that the calculated sizes match those in the output file (once strings in JSON keys are converted back to ints)
# sizes2 = dict([(int(size), count) for size, count in hn_output["sizes"].iteritems()])
# print "Sizes from output and calculation match (expect True): %s" % (sizes == sizes2)

# get command-line parameters used for the HotNet run
parameters = hn_output["parameters"]

# print out whether this was a classic or directed HotNet run
# runType = "classic" if parameters["classic"] else "directed"
# print "Was run classic or directed? (expect directed): %s" % runType 

# # load the heat scores used for the HotNet run and show them as a histogram 
heat = hnio.load_heat_json(parameters["heat_file"])
# plt.hist(heat.values(), bins=100)
# plt.show()

# get parameters used for heat file generation
heat_parameters = hn_output["heat_parameters"]

# print out the number of samples with SNVs and CNAs
genes = hnio.load_genes(heat_parameters["gene_file"])
samples = hnio.load_samples(heat_parameters["sample_file"])
snvList = hnio.load_snvs(heat_parameters["snv_file"], genes, samples)
cnaList = hnio.load_cnas(heat_parameters["cna_file"], genes, samples)
# print "There were %s samples with SNVs" % (len(snvList))
# print "There were %s samples with CNAs" % (len(cnaList))

# check whether the first two genes in the largest CC interact
edges = hnio.load_ppi_edges(parameters["edge_list_file"])
index2gene = hnio.load_index(parameters["infmat_index_file"])
gene2index_mapping = {gene: index for index, gene in index2gene.items()}
gene1index = gene2index_mapping[components[0][0]]
gene2index = gene2index_mapping[components[0][1]]
# genes_interact = (gene1index, gene2index) in edges or (gene2index, gene1index) in edges
# print "Do the first two genes in the largest CC interact (expect False)? %s" % (genes_interact)

# how about hte first gene and the third gene?
gene3index = gene2index_mapping[components[0][2]]
genes_interact = (gene1index, gene3index) in edges or (gene3index, gene1index) in edges
# print "Do the first and third genes in the largest CC interact (expect True)? %s" % (genes_interact)


query = ["network", "oncoprint", "loliplot"]

# This loads the annotation data
annotation_data = json.load(open("/gpfs/main/home/bournewu/pan12.json"))["genes"]

# This remembers the association of each sample to every gene that it contains
sample_dict = {sample: {} for sample in samples}

# This remembers the association of each gene with every sample that contains it
gene_dict = {gene: {} for gene in genes}

cancer_dict = {}

cancer_gene_file = [line.split(" ") for line in open("samples.lst","r")]
for i in range(len(cancer_gene_file)):
	littlelist=cancer_gene_file[i]
	cancer_dict[littlelist[0]]=littlelist[1].split()

# with open('hotnet_viz_data.json') as data_file:
# 	data=json.load(data_file)

def parseInactivatingGenes(mutationList = "temp/inactivating_genes.tsv"):
	mutation_list = [line.split() for line in open(mutationList)]
	mutation_dict = {}
	for line in mutation_list:
		sample_name = line[0]
		mutation_dict[sample_name] = []
		for i in range(1, len(line)):
			mutation_dict[sample_name].append(line[i])
	return mutation_dict

inactivating_genes = parseInactivatingGenes()

def parseMutationList(mutationList, category):
	for sample in mutationList:
		cur_sample = str(sample.sample)
		cur_gene = str(sample.gene)
		if (cur_sample not in gene_dict[cur_gene]):
			gene_dict[cur_gene][cur_sample] = cur_sample
		inactivating = False
		if (cur_sample in inactivating_genes):
			inactivating = True if (cur_gene in inactivating_genes[cur_sample]) else False
		if (cur_gene not in sample_dict[cur_sample]):
			sample_dict[cur_sample][cur_gene] = { 
			"gene": cur_gene,
			"sample": cur_sample,
			"cancer": cancer_dict[cur_sample],
			"inactivating": inactivating}
		sample_dict[cur_sample][cur_gene][category] = sample.mut_type

parseMutationList(snvList, "snv")
parseMutationList(cnaList, "cna")


# Inputs:	return_list: 		a list holding all the dictionaries that
# 								holds the dictionaries representing the 
#								json for each subnetwork
# 			components_list:	the list of list of genes in each subnetwork

# Output:	the return_list with the node and link fields necessary to
# 			generate the subnetworks filled in

def generateNetworkData(return_list, components_list = components):
	# Generating the nodes and links requires 
	# the links index into the respective nodes
	# This keeps track of each node's index
	index_dict = {}	
	for i in range(len(components_list)):
		nodes = []
		links = []
		index = 0
		current_subnetwork = components_list[i]

		for j in range(len(current_subnetwork)):
			gene1 = str(current_subnetwork[j])
			nodes.append({	"gene": gene1, 
							"category": i,
							"heat": heat[0][gene1]})
			index_dict[gene1] = index
			index += 1

		# Having populated the index_dict, we
		# can now generate all the links
		for j in range(len(current_subnetwork)):
			gene1 = str(current_subnetwork[j])
			for k in range(len(current_subnetwork)):
				gene2 = str(current_subnetwork[k])
				link = {"source": index_dict[gene1],
						"target": index_dict[gene2]}
				link["present"] =  True if ((gene2index_mapping[gene1], gene2index_mapping[gene2]) in edges) else False
				links.append(link)

		return_list[i]["nodes"] = nodes;
		return_list[i]["links"] = links;

	return return_list

# Inputs:	return_list: 		a list holding all the dictionaries that
# 								holds the dictionaries representing the 
#								json for each subnetwork
# 			components_list:	the list of list of genes in each subnetwork
#			sample_set:			the set of all samples
#			gene_list:			the set of all genes

# Output:	the return_list with the node and link fields necessary to
# 			generate the subnetworks filled in
def generateOncoprintData(return_list,  components_list = components):

	for i in range(len(components_list)):
		current_subnetwork = components_list[i]
		samples = []
		for j in range(len(current_subnetwork)):
			gene1 = str(current_subnetwork[j])
			cur_sample_list = gene_dict[gene1].keys()
			for k in range(len(cur_sample_list)):
				cur_sample = cur_sample_list[k]
				filtered_gene_list = []
				for l in range(len(current_subnetwork)):
					temp_cur_component = current_subnetwork[l]
					if (temp_cur_component in sample_dict[cur_sample]):
						filtered_gene_list.append(sample_dict[cur_sample][temp_cur_component])
				samples.append({	"name": cur_sample,
									"cancer": cancer_dict[cur_sample],
									"genes": filtered_gene_list
									})
		return_list[i]["samples"] = samples
	return return_list

def generateLolliplotData(return_list, components_list = components):

	for i in range(len(components_list)):
		current_subnetwork = components_list[i]
		gene_mutations = []
		for current_gene in current_subnetwork:
			transcripts = []
			if current_gene in annotation_data:
				current_annotation_data = annotation_data[current_gene]
				for transcript_name in current_annotation_data:
					current_transcript = current_annotation_data[transcript_name]
					current_mutations = current_transcript["mutations"]
					mutation_array = []
					for mutation in current_mutations:
						locus = mutation["loc"]
						types = mutation["types"]
						for mutation_type in types:
							current_mutation_type = mutation_type["name"]
							current_samples = mutation_type["samples"]
							for current_sample in current_samples:
								mutation_array.append({	"locus": locus,
														"sample": current_sample,
														"mutation": current_mutation_type ,
														"cancer": mutation_type["cancer"] })


					transcript = {	"name": transcript_name,
									"gene": current_gene,
									"length": current_transcript["length"],
									"domains": current_transcript["domains"],
									"mutations":mutation_array }
					transcripts.append(transcript)


			gene_mutations.append({	"gene": current_gene, 
									"transcripts": transcripts})
#									"mutations": list_of_mutations})
		return_list[i]["annotations"] = gene_mutations
	return return_list
#		END OF HOTNET CODE. IT'S VISUALIZATION TIME!

return_list = [{} for x in range(len(components))]

if ("network" in query):
	return_list = generateNetworkData(return_list)

if ("oncoprint" in query):
	return_list = generateOncoprintData(return_list)

if ("loliplot" in query):
	return_list = generateLolliplotData(return_list);

with open('hotnet_viz_data.json', 'w+') as outfile:
	json.dump(return_list, outfile, skipkeys=False, ensure_ascii=True, indent=1 )