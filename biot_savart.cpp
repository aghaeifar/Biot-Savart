
#include <execution>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <cmath>
#include "biot_savart.h"

#ifdef MATLAB_MEX_FILE 
#include "mex.h"
#include "matrix.h"
#endif

const float coef1 = 4.0/3.0;
const float coef2 = 1.0/3.0;

void integrate(const float *dL, 
                const float *R,
                float *b1)
{
    float r  = sqrt(R[0]*R[0] + R[1]*R[1] + R[2]*R[2]);
    float ir = 1.f / (r * r * r);

    b1[0] = ir * (dL[1] * R[2] - dL[2] * R[1]);
    b1[1] = ir * (dL[2] * R[0] - dL[0] * R[2]);
    b1[2] = ir * (dL[0] * R[1] - dL[1] * R[0]);
}

void core(const float *start,
                const float *end,
                const float *grid,
                float *b1)
{
    float start_end[3]  = {end[0] + start[0] , end[1] + start[1] , end[2] + start[2]};
    float mid[3]        = {0.5f*start_end[0] , 0.5f*start_end[1] , 0.5f*start_end[2]};
    float mid_up[3]     = {0.75f*start_end[0], 0.75f*start_end[1], 0.75f*start_end[2]};
    float mid_down[3]   = {0.25f*start_end[0], 0.25f*start_end[1], 0.25f*start_end[2]};

    float dL_full[3]    = {end[0] - start[0], end[1] - start[1], end[2] - start[2]};
    float dL_half[3]    = {0.5f*dL_full[0], 0.5f*dL_full[1], 0.5f*dL_full[2]};

    float R[3]          = {grid[0] - mid[0]     , grid[1] - mid[1]     , grid[2] - mid[2]};
    float R_up[3]       = {grid[0] - mid_up[0]  , grid[1] - mid_up[1]  , grid[2] - mid_up[2]};
    float R_down[3]     = {grid[0] - mid_down[0], grid[1] - mid_down[1], grid[2] - mid_down[2]};

    float b1_full[3]    = {0.f};
    float b1_up[3]      = {0.f};
    float b1_down[3]    = {0.f};

    integrate(dL_full, R, b1_full);
    integrate(dL_half, R_up, b1_up);
    integrate(dL_half, R_down, b1_down);

    b1[0] += coef1*(b1_up[0] + b1_down[0]) - coef2*b1_full[0];
    b1[1] += coef1*(b1_up[1] + b1_down[1]) - coef2*b1_full[1];
    b1[2] += coef1*(b1_up[2] + b1_down[2]) - coef2*b1_full[2];
}



bool calculate_b1(const uint32_t num_seg, 
                  const float *seg_start_xyz, 
                  const float *seg_end_xyz, 
                  const uint32_t num_grid, 
                  const float *grid_xyz, 
                  float *b1_xyz)
{ 
    std::fill(b1_xyz, b1_xyz + 3*num_grid, 0.f);

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
        b1_xyz[grid] *= 1e-7;

    return true;
}


#ifdef __cplusplus
extern "C" {
#endif
bool simulate(const uint32_t num_seg, 
              const float *seg_start_xyz, 
              const float *seg_end_xyz, 
              const uint32_t num_grid, 
              const float *grid_xyz, 
              float *b1_xyz)
{ 
    return calculate_b1(num_seg, seg_start_xyz, seg_end_xyz, num_grid, grid_xyz, b1_xyz);
}
#ifdef __cplusplus
} // extern "C"
#endif


#ifdef MATLAB_MEX_FILE 

_T* getPointer(const mxArray *prhs)
{
    #ifdef __SINGLE_PRECISION__
    return mxGetSingles(prhs);
    #else
    return mxGetDoubles(prhs);
    #endif
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    if (nrhs < 3)
        mexErrMsgTxt("Wrong number of inputs.");

    for (int i=1; i<nrhs; i++)
    #ifdef __SINGLE_PRECISION__
        if (mxIsDouble(prhs[i]))
            mexErrMsgTxt("all inputs must be single!");
    #else
        if (mxIsSingle(prhs[i]))
            mexErrMsgTxt("all inputs must be double!");
    #endif

    if (mxGetN(prhs[0]) != mxGetN(prhs[1]) || mxGetM(prhs[0]) != 3 || mxGetM(prhs[1]) != 3 || mxGetM(prhs[2]) != 3)
        mexErrMsgTxt("Input dimensions are wrong. Program expect following dimensions: start,3xN; end,3xN; grid,3xM");

    _T *start = getPointer(prhs[0]);
    _T *end   = getPointer(prhs[1]);
    _T *grid  = getPointer(prhs[2]);

    mwSize dims[2] = {3, mxGetN(prhs[2])};
    plhs[0] = mxCreateNumericArray(2, dims, 
    #ifdef __SINGLE_PRECISION__
            mxSINGLE_CLASS, 
    #else
            mxDOUBLE_CLASS,
    #endif
            mxREAL);
    _T *b1 = getPointer(plhs[0]);
    if (calculate_b1(mxGetN(prhs[0]), start, end, mxGetN(prhs[2]), grid, b1) == false)
        mexErrMsgTxt("Simulation failed.");
}

#endif