
#ifndef _BIOTSAVART_
#define _BIOTSAVART_


#ifdef WIN32
#ifdef __EXPORT_CLASS_BIOTSAVART__
#define DllExport __declspec(dllexport)
#else
#define DllExport __declspec(dllimport)
#endif
#else
#define DllExport
#endif


#ifdef __SINGLE_PRECISION__
typedef float _T;
#else
typedef double _T;
#endif


void integrate(const float *dL, 
                const float *R,
                float *b1);

void core(const float *start,
                const float *end,
                const float *grid,
                float *b1);

bool calculate_b1(const uint32_t num_seg, 
                  const float *seg_start_xyz, 
                  const float *seg_end_xyz, 
                  const uint32_t num_grid, 
                  const float *grid_xyz, 
                  float *b1_xyz);

#endif // _BIOTSAVART_
