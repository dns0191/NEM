#pragma once
#include <unordered_map>
#include <vector>
#include <Eigen/Dense>

class MultiGroupNode;

// 영역(region)별 그룹(group)별 단면적 저장 (unordered_map 사용)
extern std::unordered_map<int, std::vector<std::vector<double>>> crossSections;

// 노드 저장소 (1D, 2D, 3D)
extern std::vector<MultiGroupNode*> nodeGrid1D;
extern std::vector<std::vector<MultiGroupNode*>> nodeGrid2D;
extern std::vector<std::vector<std::vector<MultiGroupNode*>>> nodeGrid3D;

class MultiGroupNode {
private:
    int number_of_groups, dim, id, region;
    int* neighbor_node[2];
    double L_l, L_r;
    std::vector<Eigen::MatrixXd> DL; // Transverse Leakage
    Eigen::VectorXd node_width, flux_avg, old_flux, new_flux, SRC, SRC1, SRC2, D_c;
    Eigen::MatrixXd A, MM, mgxs;
	std::vector<Eigen::MatrixXd> M3, M4, Q, C_m, out_current;
    MultiGroupNode* l_node, * r_node;

    void makeOneDimensionalFlux(std::vector<Eigen::MatrixXd>& C);
    void updateAverageFlux(std::vector<Eigen::MatrixXd>& C);
    void updateOutgoingCurrent(const std::vector<Eigen::MatrixXd>& C);
    void updateTransverseLeakage(int direction, int group);
    double getNodeWidth(int direction) const;
    double getSurfaceFlux(int direction, bool side, int number_of_group) const;
    double getSurfaceNetCurrent(int direction, bool side, int number_of_group) const;
    double getBeta(int direction, int number_of_group) const;
    double getDiffusionCoefficient(int number_of_group) const;
    double getIncomingCurrent(int direction, bool side, int number_of_group) const;
    double getAverageTransverseLeakage(int direction, int group) const;
    MultiGroupNode* getNeighborNode(int direction, bool side) const;

public:
    MultiGroupNode(int node_id, int node_region, int group, int dimension, const Eigen::VectorXd width);
    ~MultiGroupNode();
    void getNodeInformation() const;
    void runNEM();
    void setFluxAvg(const std::vector<double>& avgFluxValues);
    bool checkConvergence(double ERROR) const;
    int getId() const { return id; }
    int getNumberOfGroups() const { return number_of_groups; }
    double getFlux(int group) const { return flux_avg[group]; }
    double getCurrent(int dimension) const { return out_current[dimension](1, 1); }
};
