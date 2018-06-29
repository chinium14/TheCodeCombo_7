/* ======================================================================== */
/* ran_int.cpp                                                              */
/*                                                                          */
/*		This subroutine returns a random integer uniformly between range1   */
/* and range2 not including range2.											*/
/*                                                                          */
/* ======================================================================== */
double ran2();

int ran_int(int range1,int range2){
	return (int)(range1 + (range2 - range1)*ran2());
}

