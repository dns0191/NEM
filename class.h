#pragma once

#include <iostream>
#include <ostream>
#include <string>
#include <vector>

class Node {
public:
    double getNodeWidth(int direction);
    double getBeta(int direction, int group);
    double getAverageTransverseLeakage(int direction, int group);
};

class MultiGroupNode {
private:
    static int number_of_groups;
    static int dim;
    static double* DL[3][3], *node_width, * flux_avg, * old_flux, *new_flux, * out_current[3][2];
    /*
	DL: Transverse Leakage
    */
    static double** A, ** M1[3], ** M2[3], ** M3[3], ** M4[3], ** D_c, **MM;
    /*
	A: Removal Cross Section
    M1: AC1 + DL1
	M2: AC2 + DL2
    M3: (6/node_width^2 + 1/10 * A)^-1 * M1
	M4: (10/node_width^2 + 1/14 * A)^-1 * M2
	D_c: Diffusion Coefficient
    MM: Q0
    */
    static double* Q[3][4];
    static double* C0[3], * C1[3], * C2[3], * C3[3], * C4[3];
    static double* B3[3][3], * B4[3][3];
    static double* SRC, * SRC1, * SRC2;
    //SRC: Node Average Flux -> s 
    static int neighbor_node[6];
    static double mgxs;

public:
    static void makeOneDimensionalFlux(double* C[3][5], double* source_avg, double* surf_source[3][2]);
    static void updateAverageFlux(double* C[3][5], const double* source_avg);
    static void updateOutgoingCurrent(double* C[3][5]);
    static void updateTransverseLeakage(int u, int group);
    static void getDimension(double* matrix);
    static double getNodeWidth(int direction);
    static double getSurfaceFlux(int u, bool side, unsigned int number_of_group);
    static double getSurfaceNetCurrent(int direction, bool side, unsigned int number_of_group);
    static double getBeta(double direction, unsigned int number_of_group);
    static double getDiffusionCoefficient(unsigned int number_of_group);
    static double getIncomingCurrent(int direction, bool side, unsigned int number_of_group);
    static Node* getNeighborNode(int direction, bool side);
    //static void Initializer(const std::string&filename);
    //static double getDiffer(double* flux1, double* flux2);
};

