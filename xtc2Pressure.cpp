// calculate the pressure tensor
// read trr or xtc files
// output Pxx(z), Pyy(z) and Pzz(z)
//
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cstring>
#include <fstream>
#include <ctime>
#include <malloc.h>
#include </usr/include/xdrfile/xdrfile_trr.h> // xdr include file 
#include </usr/include/xdrfile/xdrfile_xtc.h> // xdr include file 

    using namespace std;

double getZs(rvec *coor,matrix box,double dx,int natoms,int natomPmol);

int main(int argc, char **argv){
	
	/*variables*/
	// XTC variables
	XDRFILE *xd;		// the xtc file

	int natoms;	// number of total atoms
	int step;		// the step counter for the simulation
	float time;		// simulation time

	matrix box;		// box coordinates in a 3x3 matrix
	rvec *coor;		// atom coordinates in a 2-D matrix
	rvec *vel;		// atom coordinates in a 2-D matrix
	rvec *force;		// atom coordinates in a 2-D matrix
//	float prec;	
	float lambda;
	
	int firstframe = 500;
	if(argc >= 3) firstframe = atoi(argv[2]);
	printf("# first frame = %d \n",firstframe);

	int lastframe = 2000;			// ps
	if(argc >= 4) lastframe = atoi(argv[3]);
	printf("# last frame = %d \n",lastframe);

	int natomPmol = 1;
	if(argc >= 5) natomPmol = atoi(argv[4]);
	printf("# atoms per molecules= %d \n",natomPmol);

	double binSize = 0.1;
	if(argc >= 6) binSize = atof(argv[5]);
	printf("# binSize = %f \n",binSize);

	double temperature = 148.0;
	if(argc >= 7) temperature = atof(argv[6]);
	printf("# temperature = %f \n",temperature);

	int nCut = 1000;
	if(argc >= 8) nCut= atoi(argv[7]);
	printf("# nCut = %d \n",nCut);

	//other variables
	
	double zs;

	int nframe = lastframe+1;
	int zbin;	// divide box into zbin slabs.
	
	double dx,dy,dz,dr;	// distance between i and j
	double fij;		// force between i and j
	double vSlab;
	int nCount = 0;	// count number of frames

	double epsilon = 1.23047;
	double sigma = 0.373;
	double rCut = 2.50*sigma;
	double beta = 1000.0/8.314/temperature;

	double **pzz;
	double **pxx;
	double **pyy;
		pzz = (double **)malloc(sizeof(double *)*nframe);	
		pxx = (double **)malloc(sizeof(double *)*nframe);	
		pyy = (double **)malloc(sizeof(double *)*nframe);	

	double pzzTime;

	//char f[40];
	FILE *data;
	int i,j;

	/*read xtc files*/
	read_trr_natoms(argv[1],&natoms);

	coor = (rvec *)malloc(natoms*sizeof(rvec));
	vel = (rvec *)malloc(natoms*sizeof(rvec));
	force = (rvec *)malloc(natoms*sizeof(rvec));

	//open xtc file and loop through each frame
	xd=xdrfile_open(argv[1],"r");
	while( ! read_trr(xd, natoms, &step, &time, &lambda, box, coor, vel, force)){
		if(coor == 0){
			printf("Insufficient memory to load .trr file. \n");
			return 1;
		}
		if(step == 0){
			zbin = int(box[2][2]/binSize);
			vSlab = binSize*box[0][0]*box[1][1];
			for(i=0;i<nframe;i++){
				pzz[i] = (double *)malloc(sizeof(double)*zbin);
				pxx[i] = (double *)malloc(sizeof(double)*zbin);
				pyy[i] = (double *)malloc(sizeof(double)*zbin);
				for(j=0;j<zbin;j++){
					pzz[i][j] = 0.0;
					pxx[i][j] = 0.0;
					pyy[i][j] = 0.0;
				}
			}
			printf("# zbin = %d\n",zbin);
			data = fopen("pzz_Time.dat","w");
		}

	    if(time >= 0 && ( time <= lastframe || lastframe < firstframe ) ){
			
			/*** get zs to align pressure tensors in each frame **/
			zs = getZs(coor,box,binSize,natoms,natomPmol);		// get comZ in each configuration
			if(zs > box[2][2]){
				return 1;
			}

			/*** calculate pressure tensor ***/
			// ideal part
			for(i=0;i<natoms/natomPmol;i++){
				if(coor[i*natomPmol][2] < zs){
					pzz[nCount][ int((coor[i*natomPmol][2] - zs + box[2][2])/binSize) ] += 1.0/beta;
					pxx[nCount][ int((coor[i*natomPmol][2] - zs + box[2][2])/binSize) ] += 1.0/beta;
					pyy[nCount][ int((coor[i*natomPmol][2] - zs + box[2][2])/binSize) ] += 1.0/beta;
				}
				if(coor[i*natomPmol][2] >= zs){
					pzz[nCount][ int((coor[i*natomPmol][2] - zs )/binSize) ] += 1.0/beta;
					pxx[nCount][ int((coor[i*natomPmol][2] - zs )/binSize) ] += 1.0/beta;
					pyy[nCount][ int((coor[i*natomPmol][2] - zs )/binSize) ] += 1.0/beta;
				}
			}

			for(i=0;i<natoms/natomPmol;i++){
				for(j=i+1;j<natoms/natomPmol;j++){
					dx = coor[i*natomPmol][0] - coor[j*natomPmol][0];
					dy = coor[i*natomPmol][1] - coor[j*natomPmol][1];
					dz = coor[i*natomPmol][2] - coor[j*natomPmol][2];
					if(dx > box[0][0]/2.0) dx = dx - box[0][0];
					if(dx < -box[0][0]/2.0) dx = dx + box[0][0];

					if(dy > box[1][1]/2.0) dy = dy - box[1][1];
					if(dy < -box[1][1]/2.0) dy = dy + box[1][1];

					if(dz > box[2][2]/2.0) dz = dz - box[2][2];
					if(dz < -box[2][2]/2.0) dz = dz + box[2][2];

					dr = sqrt(dx*dx+dy*dy+dz*dz);
					if(dr < rCut){
						fij = 4.0*epsilon/dr*(-12.0 * pow(sigma/dr,12) + 6.0*pow(sigma/dr,6) );

						if(coor[i*natomPmol][2] < zs){
							pzz[nCount][ int((coor[i*natomPmol][2] - zs + box[2][2])/binSize) ] -= 1.0/vSlab*dz*dz/dr*fij/2.0;
							pxx[nCount][ int((coor[i*natomPmol][2] - zs + box[2][2])/binSize) ] -= 1.0/vSlab*dx*dx/dr*fij/2.0;
							pyy[nCount][ int((coor[i*natomPmol][2] - zs + box[2][2])/binSize) ] -= 1.0/vSlab*dy*dy/dr*fij/2.0;
						}
						if(coor[i*natomPmol][2] >= zs){ 
							pzz[nCount][ int((coor[i*natomPmol][2] - zs )/binSize) ] -= 1.0/vSlab*dz*dz/dr*fij/2.0;
							pxx[nCount][ int((coor[i*natomPmol][2] - zs )/binSize) ] -= 1.0/vSlab*dx*dx/dr*fij/2.0;
							pyy[nCount][ int((coor[i*natomPmol][2] - zs )/binSize) ] -= 1.0/vSlab*dy*dy/dr*fij/2.0;
						}

						if(coor[j*natomPmol][2] < zs){
							pzz[nCount][ int((coor[j*natomPmol][2] - zs + box[2][2])/binSize) ] -= 1.0/vSlab*dz*dz/dr*fij/2.0;
							pxx[nCount][ int((coor[j*natomPmol][2] - zs + box[2][2])/binSize) ] -= 1.0/vSlab*dx*dx/dr*fij/2.0;
							pyy[nCount][ int((coor[j*natomPmol][2] - zs + box[2][2])/binSize) ] -= 1.0/vSlab*dy*dy/dr*fij/2.0;
						}
						if(coor[j*natomPmol][2] >= zs){ 
							pzz[nCount][ int((coor[j*natomPmol][2] - zs )/binSize) ] -= 1.0/vSlab*dz*dz/dr*fij/2.0;
							pxx[nCount][ int((coor[j*natomPmol][2] - zs )/binSize) ] -= 1.0/vSlab*dx*dx/dr*fij/2.0;
							pyy[nCount][ int((coor[j*natomPmol][2] - zs )/binSize) ] -= 1.0/vSlab*dy*dy/dr*fij/2.0;
						}

					}
				}
			} // finish loop all pairs

			pzzTime = 0.0;
			for(i=0;i<zbin-1;i++) pzzTime += 0.5*(pzz[nCount][i] + pzz[nCount][i+1])*binSize;
			fprintf(data,"%f\t%e\n",time,pzzTime);
			if(step%10000 == 0) printf("%f\t%e\n",time,pzzTime);
			nCount++;
	    }
	}
	
	printf("# number of frame counted = %d\n",nCount);	

	double pzzAverage[zbin];	
	double pxxAverage[zbin];	
	double pyyAverage[zbin];	

	for(j=0;j<zbin;j++){
		for(i=0+nCut;i<nCount;i++){
			pzzAverage[j] += pzz[i][j]/double(nCount-nCut);
			pxxAverage[j] += pxx[i][j]/double(nCount-nCut);
			pyyAverage[j] += pyy[i][j]/double(nCount-nCut);
		}
	}
	double surfaceTension = 0.0;
	for(i=0;i<zbin-1;i++){
		surfaceTension += dx*(pzzAverage[i]-0.5*(pxxAverage[i]+pyyAverage[i]) + pzzAverage[i+1]-0.5*(pxxAverage[i+1]+pyyAverage[i+1]))*0.50;
	}

	data = fopen("xtc2Pressure.dat","w");
	fprintf(data,"z\tpzz\tpxx\tpyy\tst\n");
	for(i=0;i<zbin;i++){
		fprintf(data,"%f\t%e\t%e\t%e\t%e\n",i*binSize,pzzAverage[i],pxxAverage[i],pyyAverage[i],pzzAverage[i]-0.5*(pxxAverage[i]+pyyAverage[i]));
	}

	return 0;
}

double getZs(rvec *coor,matrix box,double binSize,int natoms,int natomPmol){

	int i,j;
	double zs;
	int zbin = (box[2][2]/binSize);

	double comZ = 0.0;
	for(i=0;i<natoms/natomPmol;i++) comZ += coor[i*natomPmol][2]/double(natoms/natomPmol);

	// calculate zs
	if(abs(comZ - box[2][2]/2.0) > binSize ){
		double indicator1,indicator2;
		double integral;

		// get rho'(z)
		double *rhoZtmp;
		rhoZtmp = (double *)malloc(sizeof(double)*zbin);
		for(i=0;i<zbin;i++) rhoZtmp[i] = 0.0;

		for(i=0;i<natoms/natomPmol;i++){
			rhoZtmp[int(coor[i*natomPmol][2]/binSize)] += 1.0/double(natoms/natomPmol)/binSize;
		}
	
		indicator2 = 0.0;	
		for(i=1;i< ( zbin*10.0 );i++){
			zs = i*binSize/10.0;
			integral = 0.0;
			for(j = 0 ; (j+1)*binSize < zs ; j++)  integral += (rhoZtmp[j]+rhoZtmp[j+1])/2.0*binSize;
			integral += rhoZtmp[j]*(zs-j*binSize) ;
			indicator1 = integral - zs / box[2][2] + comZ / box[2][2] - 0.50;
			if(  indicator2*indicator1 <= 0 && (zs-box[2][2]/2.0)*(comZ - box[2][2]/2.0) <= 0 && rhoZtmp[int(zs/binSize)] < 1.2/double(zbin)/binSize  ){
				zs = binSize/10.0 / (abs(indicator1)+abs(indicator2))*abs(indicator2) + (i-1)*binSize/10.0;
				return zs;
			}
			indicator2 = indicator1;
		}

		if(indicator2*indicator1 > 0 ){
			printf("# Error when remove c.o.m...%f.(%f)..\n",comZ,box[2][2]/2.0);
			FILE *data;
			data = fopen("getZS.log","w");
			for(int zi=0;zi<zbin;zi++) fprintf(data,"%f\t%f\n",zi*binSize,rhoZtmp[zi]);

			indicator2 = 0.0;	
			for(int i = 1;i< ( zbin*10.0 );i++){
				zs = i*binSize/10.0;
				integral = 0.0;
				for(j = 0 ; (j+1)*binSize < zs ; j++) integral += (rhoZtmp[j]+rhoZtmp[j+1])/2.0*binSize;
				integral += rhoZtmp[j]*(zs-j*binSize) ;
				indicator1 = integral - zs / box[2][2] + comZ / box[2][2] - 0.50;
				zs = binSize/10.0 / (abs(indicator1)+abs(indicator2))*abs(indicator2) + (i-1)*binSize/10.0;
				fprintf(data,"%f\t%f\t%f\n",zs,indicator1,indicator2);
				indicator2 = indicator1;	
			}
			return box[2][2]+1.0;
		}
	}
	else zs = 0.0;
	return zs;
}

