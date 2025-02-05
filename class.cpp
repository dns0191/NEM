#include "class.h"

void MultiGroupNode::makeOneDimensionalFlux(double* C[3][5], double* source_avg, double* surf_source[3][2])
{
	const int ng = number_of_groups;
	const int dim = getDimension();
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
			double flux_l = getSurfaceFlux(u, LEFT_SIDE, g);
			double flux_r = getSurfaceFlux(u, RIGHT_SIDE, g);
			C[u][0][g] = flux_avg[g];
			C[u][1][g] = (flux_r - flux_l) / 2.0;
			C[u][2][g] = (flux_r + flux_l) / 2.0 - flux_avg[g];
			SRC1[g] = DL[u][1][g];
			SRC2[g] = DL[u][2][g];
		}
		add_product(SRC1, SRC2, M1[u], C[u][1], ng);
		add_product(SRC2, SRC2, M2[u], C[u][2], ng);
		GaussianElimination(M3[u], C[u][3], SRC1, ng);
		GaussianElimination(M4[u], C[u][4], SRC2, ng);
	}
}

double MultiGroupNode::updateAverageFlux(double* C[3][5], double* source_avg)
{
	const int ng = number_of_groups;
	int dim = getDimension();

	for (tg = 0; tg < ng; tg++)
	{
		for (sg = 0; sg < ng; sg++)
		{
			SRC[tg] = source_avg[tg];
		}
	}

	for (int u = 0; u < dim; u++)
	{
		node_width[u] = getNodeWidth(u);
		for (g=0; g<ng; g++)
		{
			MM[g][g] += 12.0 * Q[u][0][g] / node_width[u];
		}
	}

	for (g = 0; g < ng; g++)
	{
		for (int u = 0; u < dim; u++)
		{
			double j_in_l = getIncomingCurrent(u, LEFT_SIDE, g);
			double j_in_r = getIncomingCurrent(u, RIGHT_SIDE, g);
			double Q4 = 1.0 - Q[u][2][g] - Q[u][3][g];
			SRC[g] += (2.0 * Q[u][0][g] * C[u][4][g] + Q4*(j_in_l+j_in_r))/node_width[u];
		}
	}
	GaussianElimination(MM, new_flux, SRC, ng);
}

void MultiGroupNode::updateOutgoingCurrent(double* C[3][5])
{
	const int dim = getDimension();
	for (double u = 0; u < dim; u++)
	{
		for (unsigned int g = 0; g < number_of_groups; g++)
		{
			double j_in_l = getIncomingCurrent(u, LEFT_SIDE, g);
			double j_in_r = getIncomingCurrent(u, RIGHT_SIDE, g);
			double j_out_l = Q[u][0][g] * (6 * flux_avg[g] - C[u][4][g]) + Q[u][1][g] * C[u][3][g] + Q[u][2][g] * j_in_r + Q[u][3][g] * j_in_l;
			double j_out_r = Q[u][0][g] * (6 * flux_avg[g] - C[u][4][g]) + Q[u][1][g] * C[u][3][g] + Q[u][2][g] * j_in_l + Q[u][3][g] * j_in_r;
		}
	}
}

void MultiGroupNode::updateTransverseLeakage(int u, unsigned int g)
{
	DL[u][0][g] = 0.0;
	for (int i = 1; i < dim; i++)
	{
		int v = (u + i) % dim;
		double node_width = getNodeWidth(v);
		DL[u][0][g] += (getSurfaceNetCurrent(v.RIGHT_SIDE, g) - getSurfaceNetCurrent(v.LEFT_SIDE, g)) / node_width;
	}

	l_node = getNeighborNode(u, LEFT_SIDE);
	r_node = getNeighborNode(u, RIGHT_SIDE);

	const double h_c = getNodeWidth(u);
	double beta_c = getbeta(u, g);
	double D_c = getDiffusionCoefficient(g);
	const double DL0_c = DL[u][0][g];
	double h_l = l_node->getNodeWidth(u);
	double h_r = r_node->getNodeWidth(u);
	double beta_l = l_node->getbeta(u, g);
	double beta_r = r_node->getbeta(u, g);
	double DL0_l = l_node->getAverageTransverseLeakage(u, g);
	double DL0_r = r_node->getAverageTransverseLeakage(u, g);
	L_l = (DL0_l / h_l + DL0_c / h_c) / (beta_l + beta_c);
	L_r = (DL0_c / h_c + DL0_r / h_r) / (beta_c + beta_r);

	DL[u][1][g] = D_c * (L_r - L_l) / 2.0;
	DL[u][2][g] = D_c * (L_r + L_l - 2.0*DL0_c/D_c) / 2.0;
}
