#include "class.h"
#include "Function.h"
#include <queue>
#include <fstream>
#include <cinttypes>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <iostream>
#include <Eigen/Dense>

std::unordered_map<int, std::vector<std::vector<double>>> crossSections;
std::vector<MultiGroupNode*> nodeGrid1D;
std::vector<std::vector<MultiGroupNode*>> nodeGrid2D;
std::vector<std::vector<std::vector<MultiGroupNode*>>> nodeGrid3D;

constexpr auto LEFT_SIDE = false;
constexpr auto RIGHT_SIDE = true;

void MultiGroupNode::makeOneDimensionalFlux(std::vector<Eigen::MatrixXd>& C)
{
	const int ng = number_of_groups;
	SRC1.setZero();
	SRC2.setZero();

	for (int u = 0; u < dim; u++)
	{
		Eigen::VectorXd flux_l(ng), flux_r(ng);
		for (int g = 0; g < ng; g++)
		{
			flux_l[g] = getSurfaceFlux(u, LEFT_SIDE, g);
			flux_r[g] = getSurfaceFlux(u, RIGHT_SIDE, g);
		}

		C[u].row(0) = flux_avg;
		C[u].row(1) = 0.5 * (flux_r - flux_l);
		C[u].row(2) = 0.5 * (flux_r + flux_l) - flux_avg;

		SRC1 = DL[u].row(1);
		SRC2 = DL[u].row(2);

		add_product(SRC1, C[u].row(1), ng);
		add_product(SRC2, C[u].row(2), ng);

		Eigen::VectorXd C_row_3 = C[u].row(3);
		Eigen::VectorXd C_row_4 = C[u].row(4);

		updateC(M3[u], C_row_3, SRC1, ng);
		updateC(M4[u], C_row_4, SRC2, ng);

		C[u].row(3) = C_row_3;
		C[u].row(4) = C_row_4;
	}
}



void MultiGroupNode::updateAverageFlux(std::vector<Eigen::MatrixXd>& C)
{
	const int ng = number_of_groups;
	old_flux = flux_avg;
	SRC.setZero();

	MM = A;
	for (int u = 0; u < dim; u++)
	{
		const double node_width_u = node_width[u];
		MM.diagonal().array() += 12.0 * Q[u].row(0).array() / node_width_u;
	}

	for (int g = 0; g < ng; g++)
	{
		for (int u = 0; u < dim; u++)
		{
			const double j_in_l = getIncomingCurrent(u, LEFT_SIDE, g);
			const double j_in_r = getIncomingCurrent(u, RIGHT_SIDE, g);
			const double Q4 = 1.0 - Q[u](2, g) - Q[u](3, g);
			SRC[g] += (2.0 * Q[u](0, g) * C[u](4, g) + Q4 * (j_in_l + j_in_r)) / node_width[u];
		}
	}

	updateC(MM, new_flux, SRC, ng);
	flux_avg = new_flux;
}


void MultiGroupNode::updateOutgoingCurrent(const std::vector<Eigen::MatrixXd>& C) {
	for (int u = 0; u < dim; u++) {
		for (int g = 0; g < number_of_groups; g++) {
			const double j_in_l = getIncomingCurrent(u, LEFT_SIDE, g);
			const double j_in_r = getIncomingCurrent(u, RIGHT_SIDE, g);
			const double j_out_l = Q[u](0, g) * (6 * flux_avg[g] - C[u](4, g)) + Q[u](1, g) * C[u](3, g) + Q[u](2, g) * j_in_r + Q[u](3, g) * j_in_l;
			const double j_out_r = Q[u](0, g) * (6 * flux_avg[g] - C[u](4, g)) - Q[u](1, g) * C[u](3, g) + Q[u](2, g) * j_in_l + Q[u](3, g) * j_in_r;
			out_current[u](0, g) = j_out_l;
			out_current[u](1, g) = j_out_r;
		}
	}
}



void MultiGroupNode::updateTransverseLeakage(int direction, int group)
{
    DL[direction](0, group) = 0.0;
    for (int i = 1; i < dim; i++)
    {
        const int v = (direction + i) % dim;
        DL[direction](0, group) += (getSurfaceNetCurrent(v, RIGHT_SIDE, group) - getSurfaceNetCurrent(v, LEFT_SIDE, group)) / node_width[v];
    }

    const double beta_c = getBeta(direction, group);
    const double d_c = D_c[group];
    const double DL0_c = DL[direction](0, group);
    const double h_c = node_width[direction];

    l_node = getNeighborNode(direction, LEFT_SIDE);
    r_node = getNeighborNode(direction, RIGHT_SIDE);

    const double h_l = (l_node ? l_node->node_width[direction] : h_c);
    const double h_r = (r_node ? r_node->node_width[direction] : h_c);
    const double beta_l = (l_node ? l_node->getBeta(direction, group) : beta_c);
    const double beta_r = (r_node ? r_node->getBeta(direction, group) : beta_c);
    const double DL0_l = (l_node ? l_node->getAverageTransverseLeakage(direction, group) : DL0_c);
    const double DL0_r = (r_node ? r_node->getAverageTransverseLeakage(direction, group) : DL0_c);

    // 왼쪽 및 오른쪽 누출 계수(L_l, L_r) 계산
    L_l = (DL0_l / h_l + DL0_c / h_c) / (beta_l + beta_c);
    L_r = (DL0_c / h_c + DL0_r / h_r) / (beta_c + beta_r);

    // 최종 횡단 누출량 업데이트
    DL[direction](1, group) = d_c * (L_r - L_l) / 2.0;
    DL[direction](2, group) = d_c * (L_r + L_l - 2.0 * DL0_c / d_c) / 2.0;
}


double MultiGroupNode::getNodeWidth(int direction) const
{
	return node_width[direction];
}

double MultiGroupNode::getAverageTransverseLeakage(int direction, int group) const
{
	return DL[direction](0, group);
}

double MultiGroupNode::getBeta(int direction, int number_of_group) const
{
	return getDiffusionCoefficient(number_of_group) / node_width[direction];
}

double MultiGroupNode::getDiffusionCoefficient(int number_of_group) const
{
	return D_c[number_of_group];
}

double MultiGroupNode::getIncomingCurrent(int direction, bool side, int number_of_group) const
{
	MultiGroupNode* node = getNeighborNode(direction, side);
	return (node ? node->out_current[direction](!side, number_of_group) : out_current[direction](side, number_of_group));
}


double MultiGroupNode::getSurfaceFlux(int direction, bool side, int number_of_group) const
{
	const double j_in = getIncomingCurrent(direction, side, number_of_group);
	const double j_out = out_current[direction](side, number_of_group);
	return 2 * (j_in + j_out);
}

double MultiGroupNode::getSurfaceNetCurrent(int direction, bool side, int number_of_group) const {
	const double j_in = getIncomingCurrent(direction, side, number_of_group);
	const double j_out = out_current[direction](side, number_of_group);
	return j_out - j_in;
}

MultiGroupNode* MultiGroupNode::getNeighborNode(int direction, bool side) const {
	// 1D 저장소일 경우
	if (dim == 1) {
		// neighbor_node 초기값을 -1로 설정했으므로, 유효하지 않은 경우 nullptr 반환
		const int neighbor_id = neighbor_node[side][0];
		if (neighbor_id >= 0 && neighbor_id < static_cast<int>(nodeGrid1D.size())) {
			return nodeGrid1D[neighbor_id];
		}
	}
	// 2D 저장소일 경우
	else if (dim == 2) {
		const int y_size = static_cast<int>(nodeGrid2D[0].size());
		const int x = id / y_size;   // 몫
		const int y = id % y_size;   // 나머지

		if (direction == 0) { // X 방향
			const int neighbor_x = (side == LEFT_SIDE) ? x - 1 : x + 1;
			if (neighbor_x >= 0 && neighbor_x < static_cast<int>(nodeGrid2D.size())) {
				return nodeGrid2D[neighbor_x][y];
			}
		}
		else if (direction == 1) { // Y 방향
			const int neighbor_y = (side == LEFT_SIDE) ? y - 1 : y + 1;
			if (neighbor_y >= 0 && neighbor_y < static_cast<int>(nodeGrid2D[0].size())) {
				return nodeGrid2D[x][neighbor_y];
			}
		}
	}
	return nullptr;  // 이웃이 없으면 nullptr 반환
}

void MultiGroupNode::add_product(Eigen::VectorXd& src, const Eigen::VectorXd& C, int ng)
{
	src += A * C;
}

void MultiGroupNode::updateC(const Eigen::MatrixXd& M, Eigen::VectorXd& C, const Eigen::VectorXd& src, int ng)
{
	C = M.inverse() * src;
}

void MultiGroupNode::setFluxAvg(const std::vector<double>& avgFluxValues) {
	std::copy(avgFluxValues.begin(), avgFluxValues.end(), flux_avg.data());
}

MultiGroupNode::MultiGroupNode(int node_id, int node_region, int group, int dimension, Eigen::VectorXd width)
	: number_of_groups(group), dim(dimension), id(node_id), region(node_region), node_width(width)
{
	flux_avg = Eigen::VectorXd::Zero(group);
	new_flux = Eigen::VectorXd::Zero(group);
	old_flux = Eigen::VectorXd::Zero(group);
	SRC = Eigen::VectorXd::Zero(group);
	SRC1 = Eigen::VectorXd::Zero(group);
	SRC2 = Eigen::VectorXd::Zero(group);
	D_c = Eigen::VectorXd::Zero(group);
	mgxs = Eigen::MatrixXd::Zero(group, 4);

	
	const auto& xs_data = crossSections[region];
	for (int g = 0; g < group; g++) {
		for (int i = 0; i < 4; i++) {
			mgxs(g, i) = xs_data[g][i];
		}
	}
	for (int i = 0; i < group; i++)
	{
		D_c[i] = mgxs(i,0);
	}
	DL.resize(dimension, Eigen::MatrixXd::Zero(3, group));
	A = Eigen::MatrixXd::Zero(group, group);
	MM = Eigen::Matrix2Xd::Zero(group, group);
	for (int i = 0; i < group; ++i) {
		A(i, i) = mgxs(i, 1);
	}
	double k_eff = 1.0;
	if (group > 1) {
		A(0, 0) = mgxs(0, 1) - (mgxs(0, 3) / k_eff);
		A(0, 1) = -mgxs(0, 2); 
		A(1, 0) = -mgxs(1, 2); 
		A(1, 1) = mgxs(1, 1) - (mgxs(1, 3) / k_eff);
	}

	out_current.resize(dimension, Eigen::MatrixXd::Zero(2, group));

	M3.resize(dim, Eigen::MatrixXd::Zero(group, group));
	M4.resize(dim, Eigen::MatrixXd::Zero(group, group));

	for (int u = 0; u < dim; u++) {
		for (int i = 0; i < group; i++) {
			for (int j = 0; j < group; j++) {


				M3[u](i, j) = static_cast<double>(1) / 10 * A(i, j);
				M4[u](i, j) = static_cast<double>(1) / 14 * A(i, j);

				if (i == j) {
					M3[u](i, j) += 6 / (width[u] * width[u]) * mgxs(i, 0);
					M4[u](i, j) += 10 / (width[u] * width[u]) * mgxs(i, 0);
				}
			}
		}
	}

	C_m.resize(dimension, Eigen::MatrixXd::Zero(5, group));

	Q.resize(dimension, Eigen::MatrixXd::Zero(4, group));
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < group; j++) {
			const double BETA = getBeta(i, j);
			Q[i](0, j) = BETA / (1 + 12 * BETA);
			Q[i](1, j) = BETA / (1 + 4 * BETA);
			Q[i](2, j) = 8 * BETA / ((1 + 4 * BETA) * (1 + 12 * BETA));
			Q[i](3, j) = (1 - 48 * BETA * BETA) / ((1 + 4 * BETA) * (1 + 12 * BETA));
		}
	}

	neighbor_node[0] = new int[dimension];
	neighbor_node[1] = new int[dimension];

	l_node = nullptr;
	r_node = nullptr;
	L_l = 0.0;
	L_r = 0.0;
}

MultiGroupNode::~MultiGroupNode()
{
	// 동적 메모리 해제
	delete[] neighbor_node[0];
	delete[] neighbor_node[1];
}

void MultiGroupNode::getNodeInformation() const
{
	std::ofstream debugFile("debug.out", std::ios_base::app);
	if (!debugFile.is_open()) {
		std::cerr << "Unable to open debug.out file.\n";
		return;
	}

	debugFile << "Node ID: " << id << "\n";
	debugFile << "Region: " << region << "\n";
	debugFile << "Number of Groups: " << number_of_groups << "\n";
	debugFile << "Dimension: " << dim << "\n";

	debugFile << "Node Width: ";
	for (int i = 0; i < dim; ++i) {
		debugFile << node_width[i] << " ";
	}
	debugFile << "\n";

	debugFile << "Flux Average: ";
	for (int i = 0; i < number_of_groups; ++i) {
		debugFile << flux_avg[i] << " ";
	}
	debugFile << "\n";

	debugFile << std::fixed << std::setprecision(3);
	for (int g = 0; g < number_of_groups; ++g) {
		debugFile << "Group " << g + 1 << ": ";
		for (int i = 0; i < 4; ++i) {
			debugFile << mgxs(g, i) << " ";
		}
		debugFile << "\n";
	}
	debugFile << "\n";

	debugFile << "A values:\n";
	for (int i = 0; i < number_of_groups; ++i) {
		for (int j = 0; j < number_of_groups; ++j) {
			debugFile << A(i, j) << " ";
		}
		debugFile << "\n";
	}
	debugFile << "\n";

	debugFile << "D values:\n";
	for (int i = 0; i < number_of_groups; i++) {
		debugFile << mgxs(i, 0) << " ";
	}
	debugFile << "\n\n";

	debugFile << std::scientific;

	debugFile << "M1\n";
	for (int i = 0; i < number_of_groups; ++i) {
		debugFile << SRC1(i) << " ";
	}
	debugFile << "\n";

	debugFile << "M2\n";
	for (int i = 0; i < number_of_groups; ++i) {
		debugFile << SRC2(i) << " ";
	}
	debugFile << "\n";

	debugFile << "M3\n";
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < number_of_groups; j++) {
			for (int k = 0; k < number_of_groups; k++) {
				debugFile << M3[i](j, k) << " ";
			}
			debugFile << "\n";
		}
		debugFile << "\n";
	}
	
	debugFile << "M4\n";
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < number_of_groups; j++) {
			for (int k = 0; k < number_of_groups; k++) {
				debugFile << M4[i](j, k) << " ";
			}
			debugFile << "\n";
		}
		debugFile << "\n";
	}

	debugFile << "Q values:\n";
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < 4; j++) {
			debugFile << "Q[" << i << "][" << j << "]: ";
			for (int k = 0; k < number_of_groups; k++) {
				debugFile << Q[i](j, k) << " ";
			}
			debugFile << "\n";
		}
		debugFile << "\n";
	}

	debugFile << "C_m values:\n";
	for (int j = 0; j < 5; ++j) {
		debugFile << "C_m[" << j << "]: \n";
		for (int u = 0; u < dim; ++u) {
			for (int g = 0; g < number_of_groups; ++g) {
				debugFile << C_m[u](j, g) << " ";
			}
			debugFile << "\n";
		}
		debugFile << "\n";
	}
	debugFile << "\n";

	debugFile << "DL values:\n";
	for (int i = 0; i < 3; i++) {
		debugFile << "DL[" << i << "]" << ":\n";
		for (int j = 0; j < dim; j++) {
			for (int k = 0; k < number_of_groups; k++) {
				debugFile << DL[j](i, k) << " ";
			}
			debugFile << "\n";
		}
		debugFile << "\n";
	}

	debugFile << "out_current values:\n";
	for (int i = 0; i < 2; i++) {
		debugFile << "out_current[" << i << "]" << ":\n";
		for (int j = 0; j < dim; j++) {
			for (int k = 0; k < number_of_groups; k++) {
				debugFile << out_current[j](i, k) << " ";
			}
			debugFile << "\n";
		}
		debugFile << "\n";
	}

	debugFile << "Neighbor Node IDs:\n";
	for (int direction = 0; direction < dim; ++direction) {
		const MultiGroupNode* leftNeighbor = getNeighborNode(direction, LEFT_SIDE);
		const MultiGroupNode* rightNeighbor = getNeighborNode(direction, RIGHT_SIDE);
		debugFile << "Direction " << direction << ":\n";
		debugFile << "  Left Neighbor ID: " << (leftNeighbor ? leftNeighbor->id : -1) << "\n";
		debugFile << "  Right Neighbor ID: " << (rightNeighbor ? rightNeighbor->id : -1) << "\n";
	}
	debugFile << "\n\n";

	debugFile.close();
}

void MultiGroupNode::runNEM()
{
	const int ng = number_of_groups;

	// 벡터 연산을 활용한 횡단 누출량 업데이트
	for (int u = 0; u < dim; u++)
		for (int g = 0; g < ng; g++)
			updateTransverseLeakage(u, g);

	// 플럭스 계산 및 업데이트
	makeOneDimensionalFlux(C_m);
	updateAverageFlux(C_m);
	updateOutgoingCurrent(C_m);
}

bool MultiGroupNode::checkConvergence(double ERROR) const
{
	for (int g = 0; g < number_of_groups; ++g) {
		if (std::abs(flux_avg[g] - old_flux[g]) > ERROR) {
			return false;
		}
	}
	return true;
}
