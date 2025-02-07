#pragma once

class MultiGroupNode {
private:
    int number_of_groups;
    int dim;
    double* DL[3][3], * node_width, * flux_avg, * old_flux, *new_flux, * out_current[3][2];
    /*
	DL: Transverse Leakage
    */
    double** A, ** M1[3], ** M2[3], ** M3[3], ** M4[3], ** D_c, **MM;
    /*
	A: Removal Cross Section
    M1: AC1 + DL1
	M2: AC2 + DL2
    M3: (6/node_width^2 + 1/10 * A)^-1 * M1
	M4: (10/node_width^2 + 1/14 * A)^-1 * M2
	D_c: Diffusion Coefficient
    MM: Q0
    */
    double* Q[3][4];
    double* C0[3], * C1[3], * C2[3], * C3[3], * C4[3];
    double* B3[3][3], * B4[3][3];
    double* SRC, * SRC1, * SRC2;
    //SRC: Node Average Flux -> s
    int neighbor_node[6];
    double mgxs, L_l, L_r;
    MultiGroupNode* l_node, * r_node;

public:
    void makeOneDimensionalFlux(double* C[3][5], double* source_avg, double* surf_source[3][2]);
    void updateAverageFlux(double* C[3][5], const double* source_avg);
    void updateOutgoingCurrent(double* C[3][5]);
    void updateTransverseLeakage(int u, int group);
    void getDimension();
    double getNodeWidth(int direction);
    double getSurfaceFlux(int u, bool side, int number_of_group);
    double getSurfaceNetCurrent(int direction, bool side, int number_of_group);
    double getBeta(int direction, int number_of_group);
    double getDiffusionCoefficient(int number_of_group);
    double getIncomingCurrent(int direction, bool side, int number_of_group);
    double getAverageTransverseLeakage(int direction, int group);
    MultiGroupNode* getNeighborNode(int direction, bool side);
};

