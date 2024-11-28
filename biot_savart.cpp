
// ---------------------------------------------------------
// C++ implementation of solver for magnetic field due to a finite-length straight wires carrying a 1A (using Biot-Savart law)
// Author: A. Aghaeifar (ali.aghaeifar.mri [at] gmail.com)
// Date: 2024-2-02
// reference: https://physics.stackexchange.com/questions/662024/
// ---------------------------------------------------------

#include <algorithm>
#include <iostream>
#include <numeric>
#include <cmath>
#include <vector>   
#include <chrono>

const double tooSmall = 1e-12;

inline double dot(const double *a, const double *b)
{
    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

inline void cross(const double *a, const double *b, double *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

inline double norm (const double *a)
{
    return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}


void core(const double *start,
                const double *end,
                const double *grid,
                double *b1)
{
    double end_grid[3]    = {(end[0] - grid[0])  , (end[1] - grid[1])  , (end[2] - grid[2])};
    double start_grid[3]  = {(start[0] - grid[0]), (start[1] - grid[1]), (start[2] - grid[2])};
    double end_start[3]   = {(end[0] - start[0]) , (end[1] - start[1]) , (end[2] - start[2])};
    double norm_end_start = norm(end_start);
    double norm_end_grid  = norm(end_grid);
    double norm_start_grid= norm(start_grid);
    // distance from grid to the wire segment = end_grid * cos(phi)
    double dist[3];
    
    cross(end_grid, end_start, dist);
    
    dist[0] = dist[0] / norm_end_start;
    dist[1] = dist[1] / norm_end_start;
    dist[2] = dist[2] / norm_end_start;
    double distance = norm(dist);

    double cos_theta1 = dot(end_grid, end_start) / (norm_end_grid * norm_end_start);
    double cos_theta2 = dot(start_grid, end_start) / (norm_start_grid * norm_end_start);
    // absolute value of B1
    double absB1 = distance > tooSmall ? (cos_theta1 - cos_theta2) / distance : 0.;
    // direction of B1
    double dir[3];
    cross(end_grid, end_start, dir);

    // we need normalized dir
    double norm_dir = norm(dir);
    absB1 = norm_dir>tooSmall ? absB1/norm_dir : 0.;

    b1[0] += dir[0] * absB1;
    b1[1] += dir[1] * absB1;
    b1[2] += dir[2] * absB1;
}

// ---------------------------------------------------------
// C interface
// ---------------------------------------------------------

extern "C" {
bool calculate_b0(const uint32_t num_seg,      // number of wire segments
                  const double *seg_start_xyz, // head of segment [m], 3xN matrix, column-major order {x1,y1,z1,...,xm,ym,zm}
                  const double *seg_end_xyz,   // tail of segment [m], 3xN matrix, column-major order {x1,y1,z1,...,xm,ym,zm}
                  const uint32_t num_grid,     // number of spatial points where B1 is calculated
                  const double *grid_xyz,      // spatial points [m], 3xM matrix, column-major order {x1,y1,z1,...,xm,ym,zm}
                  double *b1_xyz               // output B1 field [T], 3xM matrix, column-major order {Bx1,By1,Bz1,...,Bxm,Bym,Bzm}
                  )
{ 
    // reset b1
    std::fill(b1_xyz, b1_xyz + 3*num_grid, 0.f);
    // loop over wire segments
    auto start = std::chrono::steady_clock::now();
    for (uint32_t seg = 0; seg < num_seg; seg++)
    {
        const auto *start = seg_start_xyz + 3*seg;
        const auto *end   = seg_end_xyz   + 3*seg;
        try
        {
            #pragma omp parallel for
            for(uint32_t grid=0; grid<num_grid; grid++)
                core(start, end, grid_xyz + 3*grid, b1_xyz + 3*grid);    
        }
        catch( std::exception &ex )
        {
            std::cout << "Simulation failed." << std::endl;
            std::cout << ex.what() << std::endl;
            return false;
        } 
    }

    for (uint32_t grid = 0; grid < 3*num_grid; grid++)
        b1_xyz[grid] *= 1e-7; // convert to T
    
    std::cout << "Simulation finished in " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start).count()/1000. << " second(s)."  << std::endl;

    return true;
}

} // extern "C"

// ---------------------------------------------------------
// MATLAB interface
// B = biot_savart(start, stop, xyz);
// ---------------------------------------------------------
#ifdef MATLAB_MEX_FILE 

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 3)
        mexErrMsgTxt("Wrong number of inputs.");

    for (int i=1; i<nrhs; i++) 
        if (mxIsDouble(prhs[i]) == false)
            mexErrMsgTxt("all inputs must be double!");

    if (mxGetN(prhs[0]) != mxGetN(prhs[1]) || mxGetM(prhs[0]) != 3 || mxGetM(prhs[1]) != 3 || mxGetM(prhs[2]) != 3)
        mexErrMsgTxt("Input dimensions are wrong. Program expect following dimensions: start,3xN; end,3xN; grid,3xM");

    double *start = mxGetDoubles(prhs[0]);
    double *end   = mxGetDoubles(prhs[1]);
    double *grid  = mxGetDoubles(prhs[2]);

    mwSize dims[2] = {3, mxGetN(prhs[2])};
    plhs[0] = mxCreateNumericArray(2, dims, mxDOUBLE_CLASS, mxREAL);
    double *b1 = mxGetDoubles(plhs[0]);
    if (calculate_b0((uint32_t)mxGetN(prhs[0]), start, end, (uint32_t)mxGetN(prhs[2]), grid, b1) == false)
        mexErrMsgTxt("Simulation failed.");
}

#endif