
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


#endif // _BIOTSAVART_
