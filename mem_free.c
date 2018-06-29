/* ======================================================================== */
/* mem_free.cpp                                                             */
/*                                                                          */
/* The subroutine frees the memory needed in initialization that is not     */
/* needed later in the program.                                             */
/*                                                                          */
/* Create by Thomas A. Knotts IV, May 31, 2007.                             */
/*                                                                          */
/* ======================================================================== */

#include "defines.h"

/* ************************************************************************ */
/*                                                                          */
/* ======================== Begin Subroutine ============================== */
/*                                                                          */
/* ************************************************************************ */
void mem_free(void) 
{
    free(nbfix_prop);
}