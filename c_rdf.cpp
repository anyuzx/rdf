#include <math.h>

/* 
c_rdf.cpp

c++ function that calculate the radial distribution function in a periodic conditon
c_rdf_one: calculate the RDF of the same type of atoms
c_rdf_two: calculate the RDF of two types of atoms 

input coordinates array need to be flatten which means that [[x1,y1,z1],[x2,y2,z2]] should be [x1,y1,z1,x2,y2,z2]

Arguments:
r_array: flatten array of coordinates of all atoms need to be calculated
rdf_array: Radial Distribution Function
dim_array: dimension array. Ex. dim_array for a cubic box whose side length is 1.0 would be [1.0,1.0,1.0]
n: number of the molecules. 
natmm,natmm1,natmm2: number of the atoms in one molecule.
nHist: number of bins of RDF
rHistMax: upper limit of the radius range of RDF
v: volume of simulation box

*/

void c_rdfone(double* r_array, double* rdf_array, double* dim_array, int n, int natmm, int nHist, double rHistMax, double v){
	int i, j, k, l, ihist;
	double dx, dy, dz, r;
	double binsize = rHistMax/nHist;
	double fact = v/(2 * 3.1415927 * n * (n-1) * natmm * natmm * binsize);

	for(i=0;i<nHist;i++){
		rdf_array[i] = 0;
	}

	for(i=0;i<n-1;i++){
		for(j=i+1;j<n;j++){
			for(k=0;k<natmm;k++){
				for(l=0;l<natmm;l++){
					dx = fabs(r_array[(i*natmm+k)*3] - r_array[(j*natmm+l)*3]);
					dy = fabs(r_array[1+(i*natmm+k)*3] - r_array[1+(j*natmm+l)*3]);
					dz = fabs(r_array[2+(i*natmm+k)*3] - r_array[2+(j*natmm+l)*3]);
					dx = fmin(dx,dim_array[0] - dx);
					dy = fmin(dy,dim_array[1] - dy);
					dz = fmin(dz,dim_array[2] - dz);
					r = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0));

					if (r < rHistMax)
					{
						ihist = int(r/binsize);
						rdf_array[ihist] = rdf_array[ihist] + fact/pow((ihist+0.5)*binsize,2.0);
					}
				}
			}
		}
	}
	return;
}

void c_rdftwo(double* r1_array, double* r2_array, double* rdf_array, double* dim_array, int n, int natmm1, int natmm2, int nHist, double rHistMax, double v){
	int i, j, k, l, ihist;
	double dx, dy, dz, r;
	double binsize = rHistMax/nHist;
	double fact = v/(4 * 3.1415927 * (n * n * natmm1 * natmm2 - n * natmm1 * natmm2) * binsize);

	for(i=0;i<nHist;i++){
		rdf_array[i] = 0;
	}

	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if (j == i){continue;}
			for(k=0;k<natmm1;k++){
				for(l=0;l<natmm2;l++){
					dx = fabs(r1_array[(i*natmm1+k)*3] - r2_array[(j*natmm2+l)*3]);
					dy = fabs(r1_array[1+(i*natmm1+k)*3] - r2_array[1+(j*natmm2+l)*3]);
					dz = fabs(r1_array[2+(i*natmm1+k)*3] - r2_array[2+(j*natmm2+l)*3]);
					dx = fmin(dx,dim_array[0] - dx);
					dy = fmin(dy,dim_array[1] - dy);
					dz = fmin(dz,dim_array[2] - dz);
					r = sqrt(pow(dx,2.0) + pow(dy,2.0) + pow(dz,2.0));

					if (r < rHistMax)
					{
						ihist = int(r/binsize);
						rdf_array[ihist] = rdf_array[ihist] + fact/pow((ihist+0.5)*binsize,2.0);
					}
				}
			}
		}
	}
	return;
}