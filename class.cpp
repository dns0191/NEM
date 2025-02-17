#include "class.h"
#include "Function.h"
#include <queue>
#include <stdexcept>
#include <fstream>
#include <cinttypes>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <cstring>  

std::unordered_map<int, std::vector<std::vector<double>>> crossSections;
std::vector<MultiGroupNode*> nodeGrid1D;
std::vector<std::vector<MultiGroupNode*>> nodeGrid2D;
std::vector<std::vector<std::vector<MultiGroupNode*>>> nodeGrid3D;

constexpr auto LEFT_SIDE = false;
constexpr auto RIGHT_SIDE = true;

void MultiGroupNode::makeOneDimensionalFlux(double** C[5])
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
        SRC1 = add_product(SRC1, C[u][1], ng);
        SRC2 = add_product(SRC2, C[u][2], ng);
        GaussianElimination(M3[u], C[u][3], SRC1, ng);
        GaussianElimination(M4[u], C[u][4], SRC2, ng);
    }
}

void MultiGroupNode::updateAverageFlux(double** C[5], const double* source_avg)
{
    const int ng = number_of_groups;
    int g;
    std::memcpy(old_flux, flux_avg, number_of_groups * sizeof(double));
    for (int tg = 0; tg < ng; tg++)
    {
        for (g = 0; g < ng; g++)
        {
            MM[tg][g] = A[tg][g];
        }
        SRC[tg] = source_avg[tg];
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
    std::memcpy(flux_avg, new_flux, number_of_groups * sizeof(double));
}

void MultiGroupNode::updateOutgoingCurrent(double** C[5]) const
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

    const double h_c = getNodeWidth(direction);
    const double beta_c = getBeta(direction, group);
    const double d_c = getDiffusionCoefficient(group);
    const double DL0_c = DL[direction][0][group];

    l_node = getNeighborNode(direction, LEFT_SIDE);
    r_node = getNeighborNode(direction, RIGHT_SIDE);

    if (l_node != nullptr)
    {
        const double h_l = l_node->getNodeWidth(direction);
        const double beta_l = l_node->getBeta(direction, group);
        const double DL0_l = l_node->getAverageTransverseLeakage(direction, group);
        L_l = (DL0_l / h_l + DL0_c / h_c) / (beta_l + beta_c);
    }
    else
        L_l = 0.0;

    if (r_node != nullptr)
    {
        const double h_r = r_node->getNodeWidth(direction);
        const double beta_r = r_node->getBeta(direction, group);
        const double DL0_r = r_node->getAverageTransverseLeakage(direction, group);
        L_r = (DL0_c / h_c + DL0_r / h_r) / (beta_c + beta_r);
    }
    else
        L_r = 0.0;

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
        if (r_node == nullptr)
        {
            return out_current[direction][1][number_of_group];
        }
        else
        {
            return r_node->out_current[direction][0][number_of_group];
        }
    }
    else
    {
        if (l_node == nullptr)
        {
            return out_current[direction][0][number_of_group];
        }
        else
        {
            return l_node->out_current[direction][1][number_of_group];
        }
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
    // 3D 저장소일 경우
    else if (dim == 3) {
        std::cout << "Unimplemented" << "\n";
    }

    return nullptr;  // 이웃이 없으면 nullptr 반환
}

double* MultiGroupNode::add_product(double* src, double* C, int ng) const
{
    const double* temp = product_matrix(A, C, ng);
    for (int i = 0; i < ng; i++) {
        src[i] += temp[i];
    }
    delete[] temp;
    return src;
}

void MultiGroupNode::GaussianElimination(double** M, double*& C, double* src, int ng)
{
    // 전진 소거 단계 (부분 피벗팅 추가)
    for (int i = 0; i < ng; i++) {
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
        for (int j = i + 1; j < ng; j++) {
            const double factor = M[j][i] / M[i][i];
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

void MultiGroupNode::setFluxAvg(const std::vector<double>& avgFluxValues) {
    std::copy(avgFluxValues.begin(), avgFluxValues.end(), flux_avg);
}

MultiGroupNode::MultiGroupNode(int node_id, int node_region, int group, int dimension, double* width)
{
    id = node_id;
    region = node_region;
    number_of_groups = group;
    dim = dimension;
    // 동적 메모리 할당 및 초기화
    node_width = new double[dimension];
    std::copy(width, width + dimension, node_width);
    flux_avg = new double[group];
    new_flux = new double[group];
    old_flux = new double[group];
    SRC = new double[group];
    SRC1 = new double[group];
    SRC2 = new double[group];

    std::memset(flux_avg, 0, group * sizeof(double));
    std::memset(new_flux, 0, group * sizeof(double));
    std::memset(old_flux, 0, group * sizeof(double));
    std::memset(SRC, 0, group * sizeof(double));
    std::memset(SRC1, 0, group * sizeof(double));
    std::memset(SRC2, 0, group * sizeof(double));

    mgxs = new double* [group];
    for (int i = 0; i < group; ++i) {
        mgxs[i] = new double[4];
        std::memset(mgxs[i], 0, 4 * sizeof(double));
    }

    if (crossSections.find(region) != crossSections.end()) {
        const auto& xs_data = crossSections[region];
        for (int i = 0; i < 4; i++) {
            for (int g = 0; g < group; g++) {
                mgxs[g][i] = xs_data[g][i];
            }
        }
    }
    else {
        std::cerr << "Warning: Region " << region << " not found in crossSections!" << "\n";
    }

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
        A[i][i] = mgxs[i][1];
    }
    for (int i = 0; i < group; i++)
    {
        for (int j = 0; j < group; j++)
        {
            A[i][j] += (-1) * mgxs[i][3];
            if (i != j)
            {
                A[i][j] += (-1) * mgxs[i][2];
            }
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

    M3 = new double** [dimension];
    M4 = new double** [dimension];
    for (int i = 0; i < dimension; ++i) {
        M3[i] = new double* [group];
        M4[i] = new double* [group];
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
        D_c[i] = mgxs[i][0];
        MM[i] = new double[group];
        std::memset(MM[i], 0, group * sizeof(double));
    }

    for (int u = 0; u < dim; u++)
    {
        for (int i = 0; i < group; i++)
        {
            for (int j = 0; j < group; j++)
            {
                M3[u][i][j] += static_cast<double>(1) / 10 * A[i][j];
                M4[u][i][j] += static_cast<double>(1) / 14 * A[i][j];
                if (i == j)
                {
                    M3[u][i][j] += 6 / width[u] * D_c[j];
                    M4[u][i][j] += 10 / width[u] * D_c[j];
                }
            }
        }
    }

    C_m = new double** [dimension];
    for (int i = 0; i < dimension; ++i) {
        C_m[i] = new double* [5];
        for (int j = 0; j < 5; ++j) {
            C_m[i][j] = new double[group];
            std::memset(C_m[i][j], 0, group * sizeof(double));
        }
    }

    Q = new double** [dimension];
    for (int i = 0; i < dimension; ++i) {
        Q[i] = new double* [4];
        for (int j = 0; j < 4; ++j) {
            Q[i][j] = new double[group];
            std::memset(Q[i][j], 0, group * sizeof(double));
        }
    }

    for (int i = 0; i < dimension; i++)
    {
        for (int j = 0; j < group; ++j)
        {
            double BETA = getBeta(i, j);
            Q[i][0][j] = BETA / (1 + 12 * BETA);
            Q[i][1][j] = BETA / (1 + 4 * BETA);
            Q[i][2][j] = 8 * BETA / ((1 + 4 * BETA) * (1 + 12 * BETA));
            Q[i][3][j] = (1 - 48 * BETA * BETA) / ((1 + 4 * BETA) * (1 + 12 * BETA));
        }
    }

    neighbor_node[0] = new int[dimension];
    neighbor_node[1] = new int[dimension];
    // 1D의 경우 neighbor_node를 -1로 초기화하여 유효하지 않음을 표시
    for (int i = 0; i < dimension; ++i) {
        neighbor_node[0][i] = -1;
        neighbor_node[1][i] = -1;
    }

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
        delete[] M3[i];
        delete[] M4[i];
    }
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

    for (int i = 0; i < dim; ++i) { // C_m 배열 해제
        for (int j = 0; j < 5; ++j) {
            delete[] C_m[i][j];
        }
        delete[] C_m[i];
    }
    delete[] C_m;

    if (mgxs) {
        // 수정: 반복 범위를 number_of_groups로 변경
        for (int i = 0; i < number_of_groups; ++i) {
            if (mgxs[i]) delete[] mgxs[i];
        }
        delete[] mgxs;
    }

    delete[] neighbor_node[0];
    delete[] neighbor_node[1];
}

void MultiGroupNode::getNodeInformation() const
{
    std::cout << "Node ID: " << id << "\n";
    std::cout << "Region: " << region << "\n";
    std::cout << "Number of Groups: " << number_of_groups << "\n";
    std::cout << "Dimension: " << dim << "\n";

    std::cout << "Node Width: ";
    if (node_width) {
        for (int i = 0; i < dim; ++i) {
            std::cout << node_width[i] << " ";
        }
    }
    else {
        std::cout << "null";
    }
    std::cout << "\n";

    std::cout << "Flux Average: ";
    if (flux_avg) {
        for (int i = 0; i < number_of_groups; ++i) {
            std::cout << flux_avg[i] << " ";
        }
    }
    else {
        std::cout << "null";
    }
    std::cout << "\n";

    if (mgxs) {
        std::cout << std::fixed << std::setprecision(3);
        for (int g = 0; g < number_of_groups; ++g) {
            std::cout << "Group " << g + 1 << ": ";
            for (int i = 0; i < 4; ++i) {
                std::cout << mgxs[g][i] << " ";
            }
            std::cout << "\n";
        }
    }
    else {
        std::cout << "null";
    }
    std::cout << "\n";

    std::cout << "C_m values:\n";
    for (int j = 0; j < 5; ++j) {
        std::cout << "C_m[" << j << "]: ";
        for (int u = 0; u < dim; ++u) {
            for (int g = 0; g < number_of_groups; ++g) {
                std::cout << C_m[u][j][g] << " ";
            }
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "A values:\n";
    for (int i = 0; i < number_of_groups; ++i) {
        for (int j = 0; j < number_of_groups; ++j) {
            std::cout << A[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "D values:\n";
    for (int i = 0; i < number_of_groups; i++)
    {
        std::cout << D_c[i] << " ";
    }
    std::cout << "\n\n";

    std::cout << "M3, M4:\n";
    for (int i = 0; i < dim; i++)
    {
        for (int j = 0; j < number_of_groups; j++)
        {
            for (int k = 0; k < number_of_groups; k++)
            {
                std::cout << M3[i][j][k] << " ";
                std::cout << M4[i][j][k] << " ";
                std::cout << "\n";
            }
            std::cout << "\n";
        }
    }

    std::cout << "Neighbor Node IDs:\n";
    for (int direction = 0; direction < dim; ++direction) {
        const MultiGroupNode* leftNeighbor = getNeighborNode(direction, LEFT_SIDE);
        const MultiGroupNode* rightNeighbor = getNeighborNode(direction, RIGHT_SIDE);
        std::cout << "Direction " << direction << ":\n";
        std::cout << "  Left Neighbor ID: " << (leftNeighbor ? leftNeighbor->id : -1) << "\n";
        std::cout << "  Right Neighbor ID: " << (rightNeighbor ? rightNeighbor->id : -1) << "\n";
    }
    std::cout << "\n\n";
}

void MultiGroupNode::runNEM()
{
    for (int u = 0; u < dim; u++)
    {
        for (int g = 0; g < number_of_groups; g++)
        {
            updateTransverseLeakage(u, g);
        }
    }
    makeOneDimensionalFlux(C_m);
    updateAverageFlux(C_m, *DL[0]);
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
