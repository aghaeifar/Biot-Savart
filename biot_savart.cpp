
#include <execution>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <cmath>
#include <vector>   

// ---------------------------------------------------------
// C++ implementation of magnetic field due to a finite-length straight wire carrying a 1A using Biot-Savart law
// Author: A. Aghaeifar (ali.aghaeifar.mri [at] gmail.com)
// Date: 2024-2-02
// reference: https://physics.stackexchange.com/questions/662024/
// ---------------------------------------------------------

const float tooSmall = 1e-10;

inline float dot(const float *a, const float *b)
{
    return (a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
}

inline void cross(const float *a, const float *b, float *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}

inline float norm (const float *a)
{
    return sqrtf(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}


void core(const float *start,
                const float *end,
                const float *grid,
                float *b1)
{
    float end_grid[3]   = {(end[0] - grid[0])  , (end[1] - grid[1])  , (end[2] - grid[2])};
    float start_grid[3] = {(start[0] - grid[0]), (start[1] - grid[1]), (start[2] - grid[2])};
    float end_start[3]  = {(end[0] - start[0]) , (end[1] - start[1]) , (end[2] - start[2])};
    float norm_end_start = norm(end_start);
    float norm_end_grid  = norm(end_grid);
    float norm_start_grid= norm(start_grid);
    // distance from grid to the wire segment = end_grid * cos(phi)
    float dist[3];
    
    cross(end_grid, end_start, dist);
    
    dist[0] = dist[0] / norm_end_start;
    dist[1] = dist[1] / norm_end_start;
    dist[2] = dist[2] / norm_end_start;
    float distance = norm(dist);

    float cos_theta1 = dot(end_grid, end_start) / (norm_end_grid * norm_end_start);
    float cos_theta2 = dot(start_grid, end_start) / (norm_start_grid * norm_end_start);
    // absolute value of B1
    float absB1 = distance > tooSmall ? (cos_theta1 - cos_theta2) / distance : 0.f;
    // direction of B1
    float dir[3];
    cross(end_grid, end_start, dir);

    // we need normalized dir
    float norm_dir = norm(dir);
    absB1 = norm_dir > tooSmall? absB1/norm_dir : 0.f;

    b1[0] += dir[0] * absB1;
    b1[1] += dir[1] * absB1;
    b1[2] += dir[2] * absB1;
}



bool calculate_b1(const uint32_t num_seg,     // number of wire segments
                  const float *seg_start_xyz, // head of segment [m], 3xN matrix, column-major order {x1,y1,z1,...,xm,ym,zm}
                  const float *seg_end_xyz,   // tail of segment [m], 3xN matrix, column-major order {x1,y1,z1,...,xm,ym,zm}
                  const uint32_t num_grid,    // number of spatial points where B1 is calculated
                  const float *grid_xyz,      // spatial points [m], 3xM matrix, column-major order {x1,y1,z1,...,xm,ym,zm}
                  float *b1_xyz               // output B1 field [mT], 3xM matrix, column-major order {Bx1,By1,Bz1,...,Bxm,Bym,Bzm}
                  )
{ 
    // reset b1
    std::fill(b1_xyz, b1_xyz + 3*num_grid, 0.f);
    // loop over wire segments
    for (int seg = 0; seg < num_seg; seg++)
    {
        const float *start = seg_start_xyz + 3*seg;
        const float *end   = seg_end_xyz   + 3*seg;
        try
        {
            std::vector<uint32_t> a(num_grid);
            std::iota (a.begin(), a.end(),0);
            std::for_each (std::execution::par_unseq, std::begin(a), std::end(a), [&](uint32_t grid){
                core(start, end, grid_xyz + 3*grid, b1_xyz + 3*grid);
            });      
        }
        catch( std::exception &ex )
        {
            std::cout << "Simulation failed." << std::endl;
            std::cout << ex.what() << std::endl;
            return false;
        } 
    }

    for (uint32_t grid = 0; grid < 3*num_grid; grid++)
        b1_xyz[grid] *= 1e-4; // convert to mT

    return true;
}

// ---------------------------------------------------------
// C interface
// ---------------------------------------------------------

extern "C" {

bool simulate(const uint32_t num_seg, 
              const float *seg_start_xyz, 
              const float *seg_end_xyz, 
              const uint32_t num_grid, 
              const float *grid_xyz, 
              float *b1_xyz)
{ 
    return calculate_b1(num_seg, seg_start_xyz, seg_end_xyz, num_grid, grid_xyz, b1_xyz);
}

} // extern "C"



// ---------------------------------------------------------
// MATLAB interface
// ---------------------------------------------------------
#ifdef MATLAB_MEX_FILE 

#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs < 3)
        mexErrMsgTxt("Wrong number of inputs.");

    for (int i=1; i<nrhs; i++)
        if (mxIsSingle(prhs[i]) == false)
            mexErrMsgTxt("all inputs must be single!");


    if (mxGetN(prhs[0]) != mxGetN(prhs[1]) || mxGetM(prhs[0]) != 3 || mxGetM(prhs[1]) != 3 || mxGetM(prhs[2]) != 3)
        mexErrMsgTxt("Input dimensions are wrong. Program expect following dimensions: start,3xN; end,3xN; grid,3xM");

    float *start = mxGetSingles(prhs[0]);
    float *end   = mxGetSingles(prhs[1]);
    float *grid  = mxGetSingles(prhs[2]);

    mwSize dims[2] = {3, mxGetN(prhs[2])};
    plhs[0] = mxCreateNumericArray(2, dims, mxSINGLE_CLASS, mxREAL);
    float *b1 = mxGetSingles(plhs[0]);
    if (calculate_b1(mxGetN(prhs[0]), start, end, mxGetN(prhs[2]), grid, b1) == false)
        mexErrMsgTxt("Simulation failed.");
}

#endif