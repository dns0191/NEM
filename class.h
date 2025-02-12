#pragma once
#include<unordered_map>

class MultiGroupNode;
std::unordered_map<int, MultiGroupNode*> node_map;


class MultiGroupNode {
private:
    static int number_of_groups;
    static int dim;
    int id, region;
    double*** DL, * node_width, * flux_avg, * old_flux, *new_flux, *** out_current;
    /*
	DL: Transverse Leakage
    */
    double** A, *** M1, *** M2, *** M3, *** M4, ** D_c, **MM;
    /*
	A: Removal Cross_Section
    M1: AC1 + DL1
	M2: AC2 + DL2
    M3: (6/node_width^2 + 1/10 * A)^-1 * M1
	M4: (10/node_width^2 + 1/14 * A)^-1 * M2
	D_c: Diffusion Coefficient
    MM: Q0
    */
    double*** Q;
    double** C0, ** C1, ** C2, ** C3, ** C4;
    double* SRC, * SRC1, * SRC2;
    //SRC: Node Average Flux -> s
    int* neighbor_node[2];
    double mgxs;
	double L_l, L_r;
    MultiGroupNode* l_node, * r_node;

    void makeOneDimensionalFlux(double** C[5], double* source_avg, double* surf_source[3][2]);
    void updateAverageFlux(double** C[5], const double* source_avg);
    void updateOutgoingCurrent(double** C[5]);
    void updateTransverseLeakage(int direction, int group);
    double getNodeWidth(int direction) const;
    double getSurfaceFlux(int direction, bool side, int number_of_group) const;
    double getSurfaceNetCurrent(int direction, bool side, int number_of_group) const;
    double getBeta(int direction, int number_of_group) const;
    double getDiffusionCoefficient(int number_of_group) const;
    double getIncomingCurrent(int direction, bool side, int number_of_group) const;
    double getAverageTransverseLeakage(int direction, int group) const;
    MultiGroupNode* getNeighborNode(int direction, bool side) const;
    static void GaussianElimination(double** M, double*& C, double* src, int ng);
    void add_product(double*& src, double**& M, double* C, int ng);

public:
    MultiGroupNode(int node_id, int node_region, int group, int dimension);
    ~MultiGroupNode();
};

