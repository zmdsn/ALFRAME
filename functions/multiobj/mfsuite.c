/*
The MATLAB wrapper for the test function suite fsuite.dll or fsuite.so
*/

/* OS specific function calls for shared objects and dynamic linking */ 


#ifdef _WIN32
#include <windows.h>
#include <process.h>
typedef void (WINAPI * PPROC) (double *, double *, int ,int);

#define LIBHANDLE HANDLE
#define GetProcedure GetProcAddress
#define CloseDynalink FreeLibrary
#else
#include <dlfcn.h>
#include <pthread.h>
typedef void (*PPROC) (double *, double *, int , int);

#define LIBHANDLE void *
#define GetProcedure dlsym
#define CloseDynalink dlclose
#endif

#include "mex.h"

#ifdef __STDC__
void mexFunction (int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]) 
#else
void  mexFunction (nlhs, plhs, nrhs, prhs) 
int   nlhs, nrhs;
mxArray *plhs[];
const mxArray *prhs[];
#endif
{
	int i, j, m, nx, n_obj, size;
	double *F, *X, *tmpptr;
	char *fname;
	PPROC pfcn;
	LIBHANDLE hLibrary;
	if ((nrhs < 3) || (nlhs < 1))
    {
		mexPrintf ("usage: f = mfsuite(x,n_obj,'function_name');\n");
		mexErrMsgTxt ("example: f = mfsuite(x, 2, 'OKA2');");
    }
	m = mxGetM (prhs[0]);
	nx = mxGetN (prhs[0]);
	X = mxGetPr (prhs[0]);
	tmpptr = mxGetPr (prhs[1]);
	n_obj = (int) tmpptr[0];
	
	size = mxGetNumberOfElements (prhs[2]) + 1;
	fname = mxCalloc (size, sizeof (char));
	if (mxGetString (prhs[2], fname, size) != 0)
		mexErrMsgTxt ("Could not convert string data.");
#ifdef _WIN32
	hLibrary = LoadLibrary ("./functions/multiobj/fsuite.dll");
#else
	hLibrary = dlopen ("./fsuite.so", RTLD_NOW);
#endif
	if (NULL == hLibrary)
		mexErrMsgTxt ("could not load fsuite.dll (win32) or fsuite.so (linux) file");
	pfcn = (PPROC) GetProcedure (hLibrary, fname);
	if (NULL == pfcn)
    {
		mexPrintf ("procedure %s not found in library file!", fname);
		mexErrMsgTxt ("failed to load procedure");
    }
	plhs[0] = mxCreateDoubleMatrix (m, n_obj, mxREAL);
	F = mxGetPr (plhs[0]);
	
	
	for (i = 0; i < m; i++)
    {
		pfcn (&X[i * nx], &F[i*n_obj], nx, n_obj);
    }
	CloseDynalink (hLibrary);
	mxFree (fname);
}
