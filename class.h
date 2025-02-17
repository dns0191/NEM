#pragma once
#include <unordered_map>
#include <vector>

class MultiGroupNode;

// 영역(region)별 그룹(group)별 단면적 저장 (unordered_map 사용)
extern std::unordered_map<int, std::vector<std::vector<double>>> crossSections;

// 노드 저장소 (1D, 2D, 3D)
extern std::vector<MultiGroupNode*> nodeGrid1D;
extern std::vector<std::vector<MultiGroupNode*>> nodeGrid2D;
extern std::vector<std::vector<std::vector<MultiGroupNode*>>> nodeGrid3D;

class MultiGroupNode {
private:
    int number_of_groups;
    int dim;
    int id, region;
    double*** DL, * node_width, * flux_avg, * old_flux, * new_flux, *** out_current;
    /*
        DL: Transverse Leakage
    */
    double** A, *** M3, *** M4, * D_c, ** MM;
    /*
        A: Removal Cross_Section
        M3: (6/node_width^2 + 1/10 * A)^-1 * M1
        M4: (10/node_width^2 + 1/14 * A)^-1 * M2
        D_c: Diffusion Coefficient
        MM: Q0
    */
    double*** Q;
    double*** C_m;
    double* SRC, * SRC1, * SRC2;
    // SRC: Node Average Flux -> s
    int* neighbor_node[2];
    double** mgxs;
    double L_l, L_r;
    MultiGroupNode* l_node, * r_node;

    void makeOneDimensionalFlux(double** C[5]);
    void updateAverageFlux(double** C[5], const double* source_avg);
    void updateOutgoingCurrent(double** C[5]) const;
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
    double* add_product(double* src, double* C, int ng) const;

public:
    MultiGroupNode(int node_id, int node_region, int group, int dimension, double* width);
    ~MultiGroupNode();
    void getNodeInformation() const;
    void runNEM();
    void setFluxAvg(const std::vector<double>& avgFluxValues);
    bool checkConvergence(double ERROR) const;
    int getId() const { return id; }
    int getNumberOfGroups() const { return number_of_groups; }
    double getFlux(int group) const { return flux_avg[group]; }
    double getCurrent(int group) const { return out_current[1][0][group]; }
};
