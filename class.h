#pragma once

class MultiGroupNode
{
	void makeOneDimensionalFlux(double *C[3][5], double *source_avg, double *surf_source[3][2]);
	double updateAverageFlux(double *C[3][5], double *source_avg);
	void updateOutgoingCurrent(double* C[3][5]);
	void updateTransverseLeakage(int u, unsigned int g);

	static unsigned int number_of_groups;
	static double* DL[3][3], * flux_avg, * old_flux, * out_current[3][2];
	static double** A, ** M1[3], ** M2[3], ** M3[3], ** M4[3], * D_c;
	static double* C0[3], * C1[3], * C2[3], * C3[3], * C4[3];
	static double* B3[3][3], * B4[3][3];
	static double* SRC, * SRC1, * SRC2;

	static int neighbor_node[6];

	static double mgxs;
};