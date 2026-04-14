#pragma once
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <functional>

struct area_element_t {
    int func_number;
    int x1, y1, x2, y2;
    area_element_t() {}
    area_element_t(int _func_number, int _x1, int _y1, int _x2, int _y2) {
        func_number = _func_number; x1 = _x1; y1 = _y1; x2 = _x2; y2 = _y2;
    }
};

struct area_t {
    std::vector<double> x_lines, y_lines;
    std::vector<int> x_indices, y_indices;
    std::vector<area_element_t> elems;
};

struct parameters_t {
    double lambda, gamma, beta; 
};

struct grid_t {
    std::vector<double> x, y;
};

struct slae_t {
    std::vector<double> u1, u2, di, l4, l5, b;
    std::vector<int> ig;
};

struct functionsBC {
    std::function<double(double, double)> u_exact = [](double x, double y) { return 2.0 * x + 3.0 * y + 1.0; };
    std::function<double(double, double, parameters_t)> f_rhs = [](double x, double y, parameters_t p) { 
        return p.gamma * (2.0 * x + 3.0 * y + 1.0); 
    };
    std::vector<std::function<double(double, double)>> firstBC = { [](double x, double y) { return 2.0 * x + 3.0 * y + 1.0; } };
    std::vector<std::function<double(double, double)>> thirdBC = { [](double x, double y) { return 2.0 * x + 3.0 * y + 1.0; } };
};

void read_area(std::string filename, area_t &area);
void read_params(std::string filename, parameters_t &params);
void read_grid(std::string filename, area_t &area, grid_t &grid);

int calcSubareaNumber(area_t &area, int i, int j, int &l);
bool isDummyNode(area_t &area, int i, int j);

void createSLAE(slae_t &slae, grid_t &grid, area_t &area, parameters_t &params, functionsBC &func);
void applyBoundaryConds(area_t &area, grid_t &g, slae_t &slae, parameters_t &par, functionsBC &func);
void clearDummyNodes(slae_t &slae, area_t &area, int x_n, int y_n);

double calcNorm(const std::vector<double> &b);
double calcDiscrepancy(const slae_t &slae, const std::vector<double> &u);
void method_BlockRelaxation(slae_t &slae, std::vector<double> &u, std::vector<double> &uk, double eps, int maxIter, grid_t &grid);

// ИЗМЕНЕНИЕ ЗДЕСЬ: добавлен параметр filename
void output(area_t &area, std::vector<double> &u, grid_t &grid, functionsBC &func, std::string filename);