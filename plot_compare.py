import matplotlib.pyplot as plt
import numpy as np
import sys

rdf = {}

input_filenam = ['RDF_200K.dat','RDF_230K.dat','RDF_235.5K.dat','RDF_298K.dat']

for filenam in input_filenam:
	check = False
	with open(filenam,'r') as f:
		for i, line in enumerate(f):
			if line.split(':')[0] == 'pair':
				pair = line.split('\n')[0].split(':')[1]
				if filenam not in rdf:
					rdf[filenam] = {}
				rdf[filenam][pair] = []
				check = True
			else:
				if check == False or not line.strip():
					continue
				rdf[filenam][pair].append([float(line.split()[0]),float(line.split()[1])])


#for filenam in iter(rdf):
#	fig, ax = plt.subplots()
#	for pair in iter(rdf[filenam]):
#		rdf[filenam][pair] = np.array(rdf[filenam][pair])
#		ax.plot(rdf[filenam][pair][:,0],rdf[filenam][pair][:,1],label = filenam)
#		legend = plt.legend(loc = 'upper right')

#	plt.xlabel(r'$r(\AA)$')
#	plt.ylabel(r'$g(r_{O-O})$')
#	plt.savefig(pair + '.png')
#	plt.clf()

for pair in iter(rdf[input_filenam[0]]):
	fig, ax = plt.subplots()
	for filenam in iter(rdf):
		rdf[filenam][pair] = np.array(rdf[filenam][pair])
		ax.plot(rdf[filenam][pair][:,0],rdf[filenam][pair][:,1],label = filenam)
		legend = plt.legend(loc = 'upper right')

	plt.xlabel(r'$r(\AA)$')
	plt.ylabel(r'$g_{O-O}(r)$')
	plt.savefig(pair + '.png')
	plt.clf()