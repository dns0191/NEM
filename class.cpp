#include "class.h"
#include <queue>
#include <fstream>
#include <cinttypes>
#include <vector>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

std::unordered_map<int, std::vector<std::vector<double>>> crossSections;
std::vector<MultiGroupNode*> nodeGrid1D;
std::vector<std::vector<MultiGroupNode*>> nodeGrid2D;
std::vector<std::vector<std::vector<MultiGroupNode*>>> nodeGrid3D;

const bool LEFT_SIDE = false;
const bool RIGHT_SIDE = true;

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

		SRC1 += A * C[u].row(1).transpose();
		SRC2 += A * C[u].row(2).transpose();

		Eigen::PartialPivLU<Eigen::MatrixXd> solver_M3(M3[u]);
		Eigen::PartialPivLU<Eigen::MatrixXd> solver_M4(M4[u]);

		C[u].row(3) = solver_M3.solve(SRC1);
		C[u].row(4) = solver_M4.solve(SRC2);
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
		MM.diagonal() += 12.0 * Q[u].row(0) / node_width[u];
	}
	for (int g = 0; g < ng; g++)
	{
		for (int u = 0; u < dim; u++)
		{
			const double j_in_r = getIncomingCurrent(u, RIGHT_SIDE, g);
			const double j_in_l = getIncomingCurrent(u, LEFT_SIDE, g);
			const double Q4 = 1.0 - Q[u](2, g) - Q[u](3, g);
			SRC[g] += (2.0 * Q[u](0, g) * C[u](4, g) + Q4 * (j_in_l + j_in_r)) / node_width[u];
		}
	}

	const Eigen::PartialPivLU<Eigen::MatrixXd> solver(MM);
	flux_avg = solver.solve(SRC);
}



void MultiGroupNode::updateOutgoingCurrent(const std::vector<Eigen::MatrixXd>& C) {
	for (int u = 0; u < dim; u++) {
		for (int g = 0; g < number_of_groups; g++) {
			const double j_in_r = getIncomingCurrent(u, RIGHT_SIDE, g);
			const double j_in_l = getIncomingCurrent(u, LEFT_SIDE, g);
			const double j_out_r = Q[u](0, g) * (6 * flux_avg[g] - C[u](4, g)) - Q[u](1, g) * C[u](3, g) - Q[u](2, g) * j_in_l + Q[u](3, g) * j_in_r;
			const double j_out_l = Q[u](0, g) * (6 * flux_avg[g] - C[u](4, g)) + Q[u](1, g) * C[u](3, g) - Q[u](2, g) * j_in_r + Q[u](3, g) * j_in_l;

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
		const double J_right = getSurfaceNetCurrent(v, RIGHT_SIDE, group);
		const double J_left = getSurfaceNetCurrent(v, LEFT_SIDE, group);

		DL[direction](0, group) += (J_right - J_left) / node_width[v];
	}

	// TL 계수 업데이트
	const double beta_c = getBeta(direction, group);
	const double d_c = D_c[group];
	const double DL0_c = DL[direction](0, group);
	const double h_c = node_width[direction];

	// 이웃 노드 참조
	l_node = getNeighborNode(direction, LEFT_SIDE);
	r_node = getNeighborNode(direction, RIGHT_SIDE);

	const double h_l = (l_node ? l_node->node_width[direction] : h_c);
	const double h_r = (r_node ? r_node->node_width[direction] : h_c);
	const double beta_l = (l_node ? l_node->getBeta(direction, group) : beta_c);
	const double beta_r = (r_node ? r_node->getBeta(direction, group) : beta_c);
	const double DL0_l = (l_node ? l_node->getAverageTransverseLeakage(direction, group) : DL0_c);
	const double DL0_r = (r_node ? r_node->getAverageTransverseLeakage(direction, group) : DL0_c);

	// TL 보정
	L_l = (DL0_l / h_l + DL0_c / h_c) / (beta_l + beta_c);
	L_r = (DL0_c / h_c + DL0_r / h_r) / (beta_c + beta_r);

	DL[direction](1, group) = d_c * (L_r - L_l) / 2.0;
	DL[direction](2, group) = d_c * (L_r + L_l - 2.0 * DL0_c / d_c) / 2.0;
}

double MultiGroupNode::getNodeWidth(int direction) const
{
	return node_width[direction];
}

double MultiGroupNode::getAverageTransverseLeakage(int direction, int group) const
{
	return DL[direction](0, group); //+ DL[direction](1,group)*2*direction/node_width[direction] + DL[direction](2, group)*(6*node_width[direction]*node_width[direction]-0.5);
}

double MultiGroupNode::getBeta(int direction, int number_of_group) const
{
	return getDiffusionCoefficient(number_of_group) / node_width[direction];
}

double MultiGroupNode::getDiffusionCoefficient(int number_of_group) const
{
	return D_c[number_of_group];
}

double MultiGroupNode::getIncomingCurrent(int direction, bool side, int number_of_group) const {
	MultiGroupNode* node = getNeighborNode(direction, side);
	if (!node) {
		// 경계 조건 적용  
		const auto it = boundaryConditions.find(direction);
		if (it != boundaryConditions.end() && it->second.find(side) != it->second.end()) {
			const BoundaryCondition condition = it->second.at(side);
			if (condition == REFLECTIVE) {
				return out_current[direction](side, number_of_group);
			}
			else
				return 0.0;
		}
		else
			return 0.0;
	}
	else
		return node->out_current[direction](!side, number_of_group);
}


double MultiGroupNode::getSurfaceFlux(int direction, bool side, int number_of_group) const {
	const double j_in = getIncomingCurrent(direction, side, number_of_group);
	const double j_out = out_current[direction](side, number_of_group);
	return 2 * (j_in + j_out);
}

double MultiGroupNode::getSurfaceNetCurrent(int direction, bool side, int number_of_group) const {
	const double j_in = getIncomingCurrent(direction, side, number_of_group);
	const double j_out = out_current[direction](side, number_of_group);
	if (side == LEFT_SIDE)
		return j_in - j_out;
	else
		return j_out - j_in;
}

MultiGroupNode* MultiGroupNode::getNeighborNode(int direction, bool side) const {
	// 1D 저장소일 경우
	if (dim == 1) {
		const int neighbor_id = neighbor_node[side][0];
		if (neighbor_id >= 0 && neighbor_id < static_cast<int>(nodeGrid1D.size())) {
			return nodeGrid1D[neighbor_id];
		}
	}
	// 2D 저장소일 경우
	else if (dim == 2) {
		const int y_size = static_cast<int>(nodeGrid2D[0].size());
		const int x = id / y_size;
		const int y = id % y_size;

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
	mg_xs = Eigen::MatrixXd::Zero(group, 4);


	const auto& xs_data = crossSections[region];
	for (int g = 0; g < group; g++) {
		for (int i = 0; i < 4; i++) {
			mg_xs(g, i) = xs_data[g][i];
		}
	}
	for (int i = 0; i < group; i++)
	{
		D_c[i] = mg_xs(i, 0);
	}
	DL.resize(dimension, Eigen::MatrixXd::Zero(3, group));
	A = Eigen::MatrixXd::Zero(group, group);
	MM = Eigen::Matrix2Xd::Zero(group, group);
	for (int i = 0; i < group; ++i) {
		A(i, i) = mg_xs(i, 1);
	}

	if (group > 1) {
		constexpr double k_eff = 1.0;
		A(0, 0) = mg_xs(0, 1) - (mg_xs(0, 3) / k_eff);
		A(0, 1) = -(mg_xs(1, 3) / k_eff);
		A(1, 0) = -mg_xs(1, 2);
		A(1, 1) = mg_xs(1, 1);
	}

	out_current.resize(dimension, Eigen::MatrixXd::Zero(2, group));
	for (int u = 0; u < dimension; ++u) {
		for (int g = 0; g < group; ++g) {
			out_current[u](0, g) = 0.1;  
			out_current[u](1, g) = 0.1; 
		}
	}

	M3.resize(dim, Eigen::MatrixXd::Zero(group, group));
	M4.resize(dim, Eigen::MatrixXd::Zero(group, group));

	for (int u = 0; u < dim; u++) {
		for (int i = 0; i < group; i++) {
			for (int j = 0; j < group; j++) {


				M3[u](i, j) = static_cast<double>(1) / 10 * A(i, j);
				M4[u](i, j) = static_cast<double>(1) / 14 * A(i, j);

				if (i == j) {
					M3[u](i, j) += 6 / (width[u] * width[u]) * mg_xs(i, 0);
					M4[u](i, j) += 10 / (width[u] * width[u]) * mg_xs(i, 0);
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

std::string MultiGroupNode::getBoundaryConditionString(int direction, bool side) const {
	const auto it = boundaryConditions.find(direction);
	if (it != boundaryConditions.end() && it->second.find(side) != it->second.end()) {
		const BoundaryCondition condition = it->second.at(side);
		switch (condition) {
		case REFLECTIVE:
			return "REFLECTIVE";
		case VACUUM:
			return "VACUUM";
		}
	}
	return "NONE";
}

void MultiGroupNode::getNodeInformation() const
{
	std::ofstream debugFile("debug.txt", std::ios_base::app);
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

	debugFile << "Boundary Conditions:\n";
	for (int direction = 0; direction < dim; ++direction) {
		debugFile << "  Direction " << direction << ":\n";
		debugFile << "    Left Side: " << getBoundaryConditionString(direction, LEFT_SIDE) << "\n";
		debugFile << "    Right Side: " << getBoundaryConditionString(direction, RIGHT_SIDE) << "\n";
	}
	debugFile << "\n";

	debugFile << std::fixed << std::setprecision(3);
	for (int g = 0; g < number_of_groups; ++g) {
		debugFile << "Group " << g + 1 << ": ";
		for (int i = 0; i < 4; ++i) {
			debugFile << mg_xs(g, i) << " ";
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
		debugFile << mg_xs(i, 0) << " ";
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

	debugFile << "MM values:\n";
	for (int i = 0; i < number_of_groups; ++i) {
		for (int j = 0; j < number_of_groups; ++j) {
			debugFile << MM(i, j) << " ";
		}
		debugFile << "\n";
	}
	debugFile << MM.determinant() << "\n";
	debugFile << "\n";

	debugFile << "s values:\n";
	for (int i = 0; i < number_of_groups; i++) {
		debugFile << SRC[i] << " ";
	}
	debugFile << "\n\n";

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
	for (int i = 0; i < dim; i++) {
		debugFile << "out_current[" << i << "]" << ":\n";
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < number_of_groups; k++) {
				debugFile << out_current[i](j, k) << " ";
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

void MultiGroupNode::getNodeInformationRaw() const{
	std::ofstream debugFile("debug.txt", std::ios_base::app);
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
	debugFile << "Boundary Conditions:\n";
	for (int direction = 0; direction < dim; ++direction) {
		debugFile << "  Direction " << direction << ":\n";
		debugFile << "    Left Side: " << getBoundaryConditionString(direction, LEFT_SIDE) << "\n";
		debugFile << "    Right Side: " << getBoundaryConditionString(direction, RIGHT_SIDE) << "\n";
	}
	debugFile << "\n";
	debugFile << std::fixed << std::setprecision(3);
	for (int g = 0; g < number_of_groups; ++g) {
		debugFile << mg_xs(g) << "\n";
	}
	debugFile << "\n";
	debugFile << "A values:\n";
	debugFile << A << "\n";
	debugFile << "\n";
	debugFile << "D values:\n";
	debugFile << D_c << "\n";
	debugFile << "\n";
	debugFile << std::scientific;
	debugFile << "M1\n";
	debugFile << SRC << "\n";
	debugFile << "M2\n";
	debugFile << SRC1 << "\n";
	debugFile << "M3\n";
	for (int i = 0; i < dim; i++) {
		debugFile << M3[i] << "\n";
	}
	debugFile << "M4\n";
	for (int i = 0; i < dim; i++) {
		debugFile << M4[i] << "\n";
	}
	debugFile << "Q values:\n";
	for (int i = 0; i < dim; i++) {
		debugFile << Q[i] << "\n";
	}
	debugFile << "MM values:\n";
	debugFile << MM << "\n";
	debugFile << "\n";
	debugFile << "s values:\n";
	for (int i = 0; i < number_of_groups; i++) {
		debugFile << SRC[i] << "\n";
	}
	debugFile << "\n";
	debugFile << "C_m values:\n";
	for (int i = 0; i < dim; ++i) {
		debugFile << C_m[i] << "\n";
	}
	debugFile << "\n";
	debugFile << "DL values:\n";
	for (int i = 0; i < dim; i++) {
		debugFile << DL[i] << "\n";
	}
	debugFile << "\n";
	debugFile << "out_current values:\n";
	for (int i = 0; i < dim; i++) {
		debugFile << out_current[i] << "\n";
	}
	debugFile << "\n";
	debugFile << "Neighbor Node IDs:\n";
	for (int direction = 0; direction < dim; ++direction) {
		const MultiGroupNode* leftNeighbor = getNeighborNode(direction, LEFT_SIDE);
		const MultiGroupNode* rightNeighbor = getNeighborNode(direction, RIGHT_SIDE);
		debugFile << "Direction " << direction << ":\n";
		debugFile << "  Left Neighbor ID: " << (leftNeighbor ? leftNeighbor->id : -1) << "\n";
		debugFile << "  Right Neighbor ID: " << (rightNeighbor ? rightNeighbor->id : -1) << "\n";
	}
	debugFile << "\n\n";
}

void MultiGroupNode::runNEM()
{
	for (int u = 0; u < dim; u++)
		for (int g = 0; g < number_of_groups; g++)
			updateTransverseLeakage(u, g);

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

void MultiGroupNode::setBoundaryCondition(int direction, bool side, BoundaryCondition condition) {
	boundaryConditions[direction][side] = condition;
}

int MultiGroupNode::getDimension() const
{
	return dim;
}