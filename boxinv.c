/* ======================================================================== */
/* boxinv()                                                                 */
/* This subroutine calculates the inverse of the axes matrix and the box    */
/* volume.                                                                  */
/* INPUT	: k                                                             */
/* OUTPUT	: update axesi[k][i] and box[k].vol                             */               
/* ======================================================================== */
#ifdef PR_NPT
#include "defines.h"

void boxinv(int ibox){
	int k;
	double detbox;

	k = ibox;

	axesi[k][0] = axes[k][4] * axes[k][8] - axes[k][7] * axes[k][5];
	axesi[k][1] = axes[k][2] * axes[k][7] - axes[k][1] * axes[k][8];
	axesi[k][2] = axes[k][1] * axes[k][5] - axes[k][2] * axes[k][4];
	axesi[k][3] = axes[k][5] * axes[k][6] - axes[k][3] * axes[k][8];
	axesi[k][4] = axes[k][0] * axes[k][8] - axes[k][2] * axes[k][6];
	axesi[k][5] = axes[k][2] * axes[k][3] - axes[k][0] * axes[k][5];
	axesi[k][6] = axes[k][3] * axes[k][7] - axes[k][4] * axes[k][6];
	axesi[k][7] = axes[k][1] * axes[k][6] - axes[k][0] * axes[k][7];
	axesi[k][8] = axes[k][0] * axes[k][4] - axes[k][1] * axes[k][3];

	detbox = axesi[k][0] * axes[k][0] + axesi[k][1] * axes[k][3] + axesi[k][2] * axes[k][6];

	axesi[k][0] = axesi[k][0] / detbox;
	axesi[k][1] = axesi[k][1] / detbox;
	axesi[k][2] = axesi[k][2] / detbox;
	axesi[k][3] = axesi[k][3] / detbox;
	axesi[k][4] = axesi[k][4] / detbox;
	axesi[k][5] = axesi[k][5] / detbox;
	axesi[k][6] = axesi[k][6] / detbox;
	axesi[k][7] = axesi[k][7] / detbox;
	axesi[k][8] = axesi[k][8] / detbox;

	box[k].vol = detbox;

	
}
#endif


