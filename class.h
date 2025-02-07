#pragma once

class MultiGroupNode {
private:
    static int number_of_groups;
    static int dim;
    double** DL[3], * node_width, * flux_avg, * old_flux, *new_flux, ** out_current[2];
    /*
	DL: Transverse Leakage
    */
    double** A, ** M1, ** M2, ** M3[3], ** M4[3], ** D_c, **MM;
    /*
	A: Removal Cross Section
    M1: AC1 + DL1
	M2: AC2 + DL2
    M3: (6/node_width^2 + 1/10 * A)^-1 * M1
	M4: (10/node_width^2 + 1/14 * A)^-1 * M2
	D_c: Diffusion Coefficient
    MM: Q0
    */
    double** Q[4];
    double** C0, ** C1, ** C2, ** C3, ** C4;
    double* B3[3][3], * B4[3][3];
    double* SRC, * SRC1, * SRC2;
    //SRC: Node Average Flux -> s
    int* neighbor_node[2];
    double mgxs, L_l, L_r;
    MultiGroupNode* l_node, * r_node;

public:
    void makeOneDimensionalFlux(double** C[5], double* source_avg, double* surf_source[3][2]);
    void updateAverageFlux(double** C[5], const double* source_avg);
    void updateOutgoingCurrent(double** C[5]);
    void updateTransverseLeakage(int direction, int group);
    void getDimension();
    double getNodeWidth(int direction) const;
    double getSurfaceFlux(int direction, bool side, int number_of_group) const;
    double getSurfaceNetCurrent(int direction, bool side, int number_of_group) const;
    double getBeta(int direction, int number_of_group) const;
    double getDiffusionCoefficient(int number_of_group) const;
    double getIncomingCurrent(int direction, bool side, int number_of_group) const;
    double getAverageTransverseLeakage(int direction, int group) const;
    MultiGroupNode* getNeighborNode(int direction, bool side);
};

