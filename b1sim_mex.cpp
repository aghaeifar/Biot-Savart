/*
Matlab mex code for simulation of Biot-Savart law (just calculate Z component of the field)
Programmer: Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>

First compile the program:
	mex b1sim_mex

Usage:
B = b1sim_mex(current, fov_grids);

current  : a cell describes wiring pattern
		   current{j}.start = [1, 2, 1]; 
		   current{j}.stop  = [4, 3, 3];
fov_grids: n*3 

*/

#include "mex.h"
#include <stdio.h>
#include <math.h>
#include "matrix.h"

template <typename T> 
T **AllocateDynamicArray( int nRows, int nCols)
{
      T **dynamicArray;
      dynamicArray = new T*[nRows];
      for( int i = 0 ; i < nRows ; i++ )
		dynamicArray[i] = new T [nCols];
      return dynamicArray;
}

template <typename T>
void FreeDynamicArray(T** dArray)
{
      delete [] *dArray;
      delete [] dArray;
}

double mx_norm(const double in[], int len)
{
	double sum = 0;
	for (int i =0; i<len; i++)
		sum += in[i]*in[i];
	return sqrt(sum);
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
	if(nrhs!=2)
        mexErrMsgTxt("2 arguments have to be passed!\n");
	mxArray *current;
    int div     = mxGetN(prhs[0]); // mxGetN = Number of columns in array  -> size(current) = 1*n (n*1 is not correct)
	int fov_len = mxGetM(prhs[1]); // mxGetM = Number of rows in array     -> size(FoV) = n*3
	double *fov_grids = mxGetPr(prhs[1]);
	double *startVal, *stopVal, subStSt[3], subStFV[3], a=0, denom;
	double *b      = new double[fov_len];
	double *c      = new double[fov_len];
	double *qz     = new double[fov_len];
	double *fz     = new double[fov_len](); // () = set to zero
	plhs[0] = mxCreateDoubleMatrix ((mwSize)fov_len, (mwSize)1, mxREAL);
	double *B = mxGetPr(plhs[0]);
	
	for(int i=0; i<div; i++)
	{
		current  = mxGetCell(prhs[0], i); // get position of the first ... 
		startVal = mxGetPr(mxGetField(current,0,"start"));//startVal[0],startVal[1],startVal[2]
		stopVal  = mxGetPr(mxGetField(current,0,"stop"));
        
		for (int j=0; j<3; j++)
			subStSt[j] = stopVal[j]-startVal[j];		
		
        a = 0; // 12.10.2017
		for (int ii =0; ii<3; ii++)
			a += subStSt[ii]*subStSt[ii];
        //mexPrintf("a=%f\n", a);        
		if(a != 0)
		{
			for(int j=0;j<fov_len;j++)
			{
				subStFV[0] = startVal[0]-fov_grids[j+0*fov_len];
				subStFV[1] = startVal[1]-fov_grids[j+1*fov_len];
				subStFV[2] = startVal[2]-fov_grids[j+2*fov_len];
				
				b[j] = 2*(subStSt[0]*subStFV[0]+subStSt[1]*subStFV[1]+subStSt[2]*subStFV[2]);	
				c[j] = subStFV[0]*subStFV[0]+subStFV[1]*subStFV[1]+subStFV[2]*subStFV[2];// //using the pow makes the function very slow --> pow(subStFV[0],2)+ pow(subStFV[1],2)+ pow(subStFV[2],2);	
				qz[j]= subStSt[0]*subStFV[1] - subStSt[1]*subStFV[0];
				
				//  integral_exact
				denom = 1/(4*a*c[j]-b[j]*b[j]); //change division by multiplication for the next line --> 4*a*c[j]-b[j]*b[j];
				fz[j] += 2*qz[j]*denom*( (2*a+b[j])/sqrt(a+b[j]+c[j]) - b[j]/sqrt(c[j]) ); //change division by multiplication --> qz[j]*(2*(2*a+b[j])/denom/sqrt(a+b[j]+c[j]) - 2*b[j]/denom/sqrt(c[j]));
				B[j] = fz[j] * 1e-4;
			}			
		}
	}
	delete[] b;
	delete[] c;
	delete[] qz; 
	delete[] fz;
//  mexPrintf("%d %d\n", dims[0], dims[1]);
}
