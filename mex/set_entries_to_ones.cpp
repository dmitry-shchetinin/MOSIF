/*=================================================================
* this routine creates the copy of input atrix in sparse format and 
* replaces all nonzero elements with ones
*=================================================================*/


#include "mex.h"


void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )
{
      
    double *A;
    mwIndex i;
    mwSize n;
    
    // create deep copy of input to output
    plhs[0] = mxDuplicateArray(prhs[0]);
    
    // get pointer to vector of nonzero elements in the output matrix
    A = mxGetDoubles(plhs[0]);
    
    //get the number of structurally nonzero entries in the matrix
    n = mxGetNzmax(plhs[0]);
    
    //replace nonzero elements with ones
    for ( i = 0; i < n; i++ ) 
        A[i] = 1;
}