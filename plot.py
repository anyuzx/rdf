import matplotlib.pyplot as plt
import numpy as np

rdf = {}

with open('RDF.dat','r') as f:
	for i,line in enumerate(f):
		if line.split(':')[0] == 'pair':
			pair = line.split(':')[1]
			rdf[pair] = []
		else:
			if not line.strip():
				continue
			rdf[pair].append([float(line.split()[0]),float(line.split()[1])])

for pair in iter(rdf):
	rdf[pair] = np.array(rdf[pair])
	plt.plot(rdf[pair][:,0],rdf[pair][:,1])
	plt.xlabel(r'$r(\AA)$')
	plt.ylabel(r'$g(r)$')
	plt.savefig(pair + '.png')
	plt.clf()