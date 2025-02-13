#include "class.h"
#include "Function.h"
#include <queue>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <vector>

constexpr auto LEFT_SIDE = false;
constexpr auto RIGHT_SIDE = true;

void MultiGroupNode::makeOneDimensionalFlux(double** C[5], double* source_avg, double* surf_source[3][2])
{
	const int ng = number_of_groups;
	int g;
	for (g = 0; g < ng; g++)
	{
		SRC1[g] = 0.0;
		SRC2[g] = 0.0;
	}

	for (int u = 0; u < dim; u++)
	{
		for (g = 0; g < ng; g++)
		{
			const double flux_l = getSurfaceFlux(u, LEFT_SIDE, g);
			const double flux_r = getSurfaceFlux(u, RIGHT_SIDE, g);
			C[u][0][g] = flux_avg[g];
			C[u][1][g] = (flux_r - flux_l) / 2.0;
			C[u][2][g] = (flux_r + flux_l) / 2.0 - flux_avg[g];
			SRC1[g] = DL[u][1][g];
			SRC2[g] = DL[u][2][g];
		}
		add_product(SRC1, M1[u], C[u][1], ng);
		add_product(SRC2, M2[u], C[u][2], ng);
		GaussianElimination(M3[u], C[u][3], SRC1, ng);
		GaussianElimination(M4[u], C[u][4], SRC2, ng);
	}
}

void MultiGroupNode::updateAverageFlux(double** C[5], const double* source_avg)
{
	const int ng = number_of_groups;
	int g;

	for (int tg = 0; tg < ng; tg++)
	{
		for (g = 0; g < ng; g++)
		{
			MM[tg][g] = A[tg][g];
			SRC[tg] = source_avg[tg];
		}
	}

	for (int u = 0; u < dim; u++)
	{
		node_width[u] = getNodeWidth(u);
		for (g = 0; g < ng; g++)
		{
			MM[g][g] += 12.0 * Q[u][0][g] / node_width[u];
		}
	}

	for (g = 0; g < number_of_groups; g++)
	{
		for (int u = 0; u < dim; u++)
		{
			const double j_in_l = getIncomingCurrent(u, LEFT_SIDE, g);
			const double j_in_r = getIncomingCurrent(u, RIGHT_SIDE, g);
			const double Q4 = 1.0 - Q[u][2][g] - Q[u][3][g];
			SRC[g] += (2.0 * Q[u][0][g] * C[u][4][g] + Q4 * (j_in_l + j_in_r)) / node_width[u];
		}
	}
	GaussianElimination(MM, new_flux, SRC, ng);
}

void MultiGroupNode::updateOutgoingCurrent(double** C[5])
{
	for (int u = 0; u < dim; u++)
	{
		for (int g = 0; g < number_of_groups; g++)
		{
			const double j_in_l = getIncomingCurrent(u, LEFT_SIDE, g);
			const double j_in_r = getIncomingCurrent(u, RIGHT_SIDE, g);
			out_current[u][0][g] = Q[u][0][g] * (6 * flux_avg[g] - C[u][4][g]) + Q[u][1][g] * C[u][3][g] + Q[u][2][g] * j_in_r + Q[u][3][g] * j_in_l;
			out_current[u][1][g] = Q[u][0][g] * (6 * flux_avg[g] - C[u][4][g]) + Q[u][1][g] * C[u][3][g] + Q[u][2][g] * j_in_l + Q[u][3][g] * j_in_r;
		}
	}
}

void MultiGroupNode::updateTransverseLeakage(int direction, int group)
{
	DL[direction][0][group] = 0.0;
	for (int i = 1; i < dim; i++)
	{
		const int v = (direction + i) % dim;
		DL[direction][0][group] += (getSurfaceNetCurrent(v, RIGHT_SIDE, group) - getSurfaceNetCurrent(v, LEFT_SIDE, group)) / getNodeWidth(v);
	}

	l_node = getNeighborNode(direction, LEFT_SIDE);
	r_node = getNeighborNode(direction, RIGHT_SIDE);

	const double h_c = getNodeWidth(direction);
	const double beta_c = getBeta(direction, group);
	const double d_c = getDiffusionCoefficient(group);
	const double DL0_c = DL[direction][0][group];
	const double h_l = l_node->getNodeWidth(direction);
	const double h_r = r_node->getNodeWidth(direction);
	const double beta_l = l_node->getBeta(direction, group);
	const double beta_r = r_node->getBeta(direction, group);
	const double DL0_l = l_node->getAverageTransverseLeakage(direction, group);
	const double DL0_r = r_node->getAverageTransverseLeakage(direction, group);
	L_l = (DL0_l / h_l + DL0_c / h_c) / (beta_l + beta_c);
	L_r = (DL0_c / h_c + DL0_r / h_r) / (beta_c + beta_r);

	DL[direction][1][group] = d_c * (L_r - L_l) / 2.0;
	DL[direction][2][group] = d_c * (L_r + L_l - 2.0 * DL0_c / d_c) / 2.0;
}

double MultiGroupNode::getNodeWidth(int direction) const
{
	return node_width[direction];
}

double MultiGroupNode::getAverageTransverseLeakage(int direction, int group) const
{
	return DL[direction][0][group];
}

double MultiGroupNode::getBeta(int direction, int number_of_group) const
{
	return getDiffusionCoefficient(number_of_group) / getNodeWidth(direction);
}

double MultiGroupNode::getDiffusionCoefficient(int number_of_group) const
{
	return D_c[number_of_group];
}

double MultiGroupNode::getIncomingCurrent(int direction, bool side, int number_of_group) const
{
	if (side == RIGHT_SIDE)
	{
		return r_node->out_current[direction][0][number_of_group];
	}
	else
	{
		return l_node->out_current[direction][1][number_of_group];
	}
}

double MultiGroupNode::getSurfaceFlux(int direction, bool side, int number_of_group) const
{
	
	if (side == RIGHT_SIDE)
	{
		return 2 * getSurfaceNetCurrent(direction, side, number_of_group);
	}
	else
	{
		return 2 * getSurfaceNetCurrent(direction, !side, number_of_group);
	}
}

double MultiGroupNode::getSurfaceNetCurrent(int direction, bool side, int number_of_group) const 
{
	
	
	if (side == RIGHT_SIDE)
	{
		const double j_in_r = getIncomingCurrent(direction, RIGHT_SIDE, number_of_group);
		const double j_out_r = out_current[direction][1][number_of_group];
		return j_in_r + j_out_r;
	}
	else
	{
		const double j_in_l = getIncomingCurrent(direction, LEFT_SIDE, number_of_group);
		const double j_out_l = out_current[direction][0][number_of_group];
		return j_in_l + j_out_l;
	}
}


MultiGroupNode* MultiGroupNode::getNeighborNode(int direction, bool side) const
{
	const int neighbor_id = neighbor_node[side][direction];

	if (neighbor_id != 0)
	{
		const auto it = node_map.find(neighbor_id);
		if(it != node_map.end())
		{
			return it->second;
		}
		else
		{
			return nullptr;
		}
	}
	else
	{
		return nullptr;
	}
}

void MultiGroupNode::add_product(double* src, double*& M, double* C, int ng)
{
	M = product_matrix(this->A, C, ng);
	for (int i=0; i<ng; i++)
	{
		src[i] += M[i];
	}
}

void MultiGroupNode::GaussianElimination(double** M, double*& C, double* src, int ng)
{
	// 전진 소거 단계 (부분 피벗팅 추가)
	for (int i = 0; i < ng; i++) {
		// 피벗 선택 (가장 큰 값으로 교환)
		int max_row = i;
		for (int j = i + 1; j < ng; j++) {
			if (fabs(M[j][i]) > fabs(M[max_row][i])) {
				max_row = j;
			}
		}
		if (max_row != i) {
			std::swap(M[i], M[max_row]);
			std::swap(src[i], src[max_row]);
		}

		// 피벗이 0이면 해가 없거나 무수히 많음
		if (M[i][i] == 0.0) {
			throw std::runtime_error("Singular matrix: No unique solution exists.");
		}

		// 소거 작업
		for (int j = i + 1; j < ng; j++) {
			double const factor = M[j][i] / M[i][i];
			for (int k = i; k < ng; k++) {
				M[j][k] -= factor * M[i][k];
			}
			src[j] -= factor * src[i];
		}
	}

	// 후진 대입 단계
	for (int i = ng - 1; i >= 0; i--) {
		C[i] = src[i];
		for (int j = i + 1; j < ng; j++) {
			C[i] -= M[i][j] * C[j];
		}
		C[i] /= M[i][i];
	}
}

MultiGroupNode::MultiGroupNode(int node_id, int node_region, int group, int dimension)
{
	id = node_id;
	region = node_region;
	// 동적 메모리 할당 및 초기화
	node_width = new double[dimension];
	flux_avg = new double[group];
	old_flux = new double[group];
	new_flux = new double[group];
	SRC = new double[group];
	SRC1 = new double[group];
	SRC2 = new double[group];

	DL = new double** [dimension];
	for (int i = 0; i < dimension; ++i) {
		DL[i] = new double* [3];
		for (int j = 0; j < 3; ++j) {
			DL[i][j] = new double[group];
			std::memset(DL[i][j], 0, group * sizeof(double));
		}
	}

	A = new double* [group];
	for (int i = 0; i < group; ++i) {
		A[i] = new double[group];
		std::memset(A[i], 0, group * sizeof(double));
	}

	Q = new double** [dimension];
	for (int i = 0; i < dimension; ++i) {
		Q[i] = new double* [4];
		for (int j = 0; j < 4; ++j) {
			Q[i][j] = new double[group];
			std::memset(Q[i][j], 0, group * sizeof(double));
		}
	}

	out_current = new double** [dimension];
	for (int i = 0; i < dimension; ++i) {
		out_current[i] = new double* [2];
		for (int j = 0; j < 2; ++j) {
			out_current[i][j] = new double[group];
			std::memset(out_current[i][j], 0, group * sizeof(double));
		}
	}

	M1 = new double* [dimension];
	M2 = new double* [dimension];
	M3 = new double** [dimension];
	M4 = new double** [dimension];

	for (int i = 0; i < dimension; ++i) {
		M1[i] = new double[group];
		M2[i] = new double[group];
		M3[i] = new double* [group];
		M4[i] = new double* [group];

		std::memset(M1[i], 0, group * sizeof(double));
		std::memset(M2[i], 0, group * sizeof(double));

		for (int j = 0; j < group; ++j) {
			M3[i][j] = new double[group];
			M4[i][j] = new double[group];
			std::memset(M3[i][j], 0, group * sizeof(double));
			std::memset(M4[i][j], 0, group * sizeof(double));
		}
	}

	D_c = new double[group];
	MM = new double* [group];
	for (int i = 0; i < group; ++i) {
		MM[i] = new double[group];
		std::memset(MM[i], 0, group * sizeof(double));
	}

	C0 = new double* [group];
	C1 = new double* [group];
	C2 = new double* [group];
	C3 = new double* [group];
	C4 = new double* [group];

	for (int i = 0; i < group; ++i) {
		C0[i] = new double[group];
		C1[i] = new double[group];
		C2[i] = new double[group];
		C3[i] = new double[group];
		C4[i] = new double[group];
		std::memset(C0[i], 0, group * sizeof(double));
		std::memset(C1[i], 0, group * sizeof(double));
		std::memset(C2[i], 0, group * sizeof(double));
		std::memset(C3[i], 0, group * sizeof(double));
		std::memset(C4[i], 0, group * sizeof(double));
	}

	mgxs = new double* [4];
	for (int i = 0; i < 4; ++i) {
		mgxs[i] = new double[group];
		std::memset(mgxs[i], 0, group * sizeof(double));
	}

	neighbor_node[0] = new int[dimension];
	neighbor_node[1] = new int[dimension];
	std::memset(neighbor_node[0], 0, dimension * sizeof(int));
	std::memset(neighbor_node[1], 0, dimension * sizeof(int));

	l_node = nullptr;
	r_node = nullptr;
	L_l = 0.0;
	L_r = 0.0;
}

MultiGroupNode::~MultiGroupNode()
{
	// 동적 메모리 해제
	delete[] node_width;
	delete[] flux_avg;
	delete[] old_flux;
	delete[] new_flux;
	delete[] SRC;
	delete[] SRC1;
	delete[] SRC2;

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < 3; ++j) {
			delete[] DL[i][j];
		}
		delete[] DL[i];
	}
	delete[] DL;

	for (int i = 0; i < number_of_groups; ++i) {
		delete[] A[i];
	}
	delete[] A;

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < number_of_groups; ++j) {
			delete[] M3[i][j];
			delete[] M4[i][j];
		}
		delete[] M1[i];
		delete[] M2[i];
		delete[] M3[i];
		delete[] M4[i];
	}

	delete[] M1;
	delete[] M2;
	delete[] M3;
	delete[] M4;

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < 4; ++j) {
			delete[] Q[i][j];
		}
		delete[] Q[i];
	}
	delete[] Q;

	for (int i = 0; i < dim; ++i) {
		for (int j = 0; j < 2; ++j) {
			delete[] out_current[i][j];
		}
		delete[] out_current[i];
	}
	delete[] out_current;

	for (int i = 0; i < number_of_groups; ++i) {
		delete[] MM[i];
	}
	delete[] D_c;
	delete[] MM;

	for (int i = 0; i < number_of_groups; ++i) {
		delete[] C0[i];
		delete[] C1[i];
		delete[] C2[i];
		delete[] C3[i];
		delete[] C4[i];
	}
	delete[] C0;
	delete[] C1;
	delete[] C2;
	delete[] C3;
	delete[] C4;

	for (int i = 0; i < 4; ++i) {
		delete[] mgxs[i];
	}
	delete[] mgxs;

	delete[] neighbor_node[0];
	delete[] neighbor_node[1];
}
