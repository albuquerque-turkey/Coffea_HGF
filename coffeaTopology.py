import sys,dendropy

def reroot_on_outgroup(tree,outgroupLeaf):
	outgroup_node = tree.find_node_with_taxon_label(outgroupLeaf)
	tree.reroot_at_edge(outgroup_node.edge, update_bipartitions=True)
	return tree


def testTopology(treeFile):
	tree=dendropy.Tree.get(path=treeFile,schema='newick',preserve_underscores=True)
	fileSplit = treeFile.split('.')
	fh = fileSplit[1]
	totalLeaves = []
	for taxon in tree.taxon_namespace:
		totalLeaves.append(taxon.label)
		if 'Gj' in taxon.label:
			outgroupLeaf = taxon.label
		elif 'Ccan' in taxon.label:
			patDipLeaf = taxon.label
		elif 'Ceug' in taxon.label:
			matDipLeaf = taxon.label
		elif 'Cara_C' in taxon.label:
			patTetLeaf = taxon.label
		elif 'Cara_E' in taxon.label:
			matTetLeaf = taxon.label
	rootedTree = reroot_on_outgroup(tree,outgroupLeaf)
	treeDict = {}
	sys.stdout.write(fh) 
	for node in rootedTree:
		subTree = node.leaf_nodes()
		currLeaves = []
		for leaf in subTree:
			if leaf.taxon.label == outgroupLeaf:
				currLeaves.append('G_jasminoides')
			elif leaf.taxon.label == patDipLeaf:
				currLeaves.append('C_canephora')
			elif leaf.taxon.label == matDipLeaf:
				currLeaves.append('C_eugenioides')
			elif leaf.taxon.label == patTetLeaf:
				currLeaves.append('C_arabica_C')
			elif leaf.taxon.label == matTetLeaf:
				currLeaves.append('C_arabica_E')
		currLeaves.sort()
		sys.stdout.write('\t' + str(currLeaves) + '\t' + str(node.label))
	sys.stdout.write('\n')


fofn = open(sys.argv[1],'r')
for line in fofn:
	testTopology(line[0:-1])

fofn.close()