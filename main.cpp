#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>
#include "fn.h"

void write_grid(int N) {
    std::ofstream out("in/grid.txt");
    if (!out.is_open()) {
        std::cerr << "Ошибка: Не удалось создать in/grid.txt!" << std::endl;
        exit(1);
    }
    out << N << " 1 " << N << " 1 " << N << " 1\n";
    out << (int)(N * 1.5) << " 1 " << N << " 1\n";
    out.close();
}

double run_test(functionsBC func, area_t &area, parameters_t &params, int N, double eps, int maxIter, std::string out_filename) {
    write_grid(N);
    grid_t grid;
    read_grid("in/grid.txt", area, grid);

    slae_t slae;
    std::vector<double> u, uk;

    createSLAE(slae, grid, area, params, func);
    applyBoundaryConds(area, grid, slae, params, func);
    clearDummyNodes(slae, area, grid.x.size(), grid.y.size());

    method_BlockRelaxation(slae, u, uk, eps, maxIter, grid);
    output(area, u, grid, func, out_filename);

    double max_err = 0;
    int n_x = grid.x.size(), n_y = grid.y.size();
    for (int j = 0; j < n_y; j++) {
        for (int i = 0; i < n_x; i++) {
            if (!isDummyNode(area, i, j)) {
                int k = j * n_x + i;
                double exact = func.u_exact(grid.x[i], grid.y[j]);
                max_err = std::max(max_err, std::abs(u[k] - exact));
            }
        }
    }
    return max_err;
}

// Новая функция: прогоняет тесты на N=10, 20, 40 и считает порядок сходимости
void run_convergence_study(std::string name, functionsBC func, area_t &area, parameters_t &params, double eps, int maxIter) {
    std::vector<int> Ns = {10, 20, 40};
    std::vector<double> errs(Ns.size());
    
    std::cout << std::left << std::setw(12) << name << " | ";
    
    for (size_t i = 0; i < Ns.size(); ++i) {
        std::string filename = "out/out_" + name + "_" + std::to_string(Ns[i]) + ".txt";
        errs[i] = run_test(func, area, params, Ns[i], eps, maxIter, filename);
        
        if (i == 0) {
            std::cout << "E_" << Ns[i] << " = " << std::scientific << std::setprecision(3) << errs[i];
        } else {
            double p = log2(errs[i - 1] / errs[i]);
            std::cout << " | p_" << Ns[i] << " = " << std::fixed << std::setprecision(3) << p;
        }
    }
    std::cout << "\n";
}

int main() {
    const double eps = 1e-13;
    const int maxIter = 50000;

    area_t area;
    parameters_t params;
    
    read_area("in/area.txt", area);
    read_params("in/params.txt", params);
    
    params.lambda = 1.0; params.gamma = 1.0; params.beta = 2.0;

    if (area.x_lines.size() < 4 || area.y_lines.size() < 3) {
        std::cerr << "КРИТИЧЕСКАЯ ОШИБКА: Файл in/area.txt не прочитался или пуст!\n";
        return 1; 
    }
    
    // --- Инициализация функций (Линейная, квадратичная и т.д.) ---
    functionsBC f_lin;
    f_lin.u_exact = [](double x, double y) { return 2.0 * x + 3.0 * y + 1.0; };
    f_lin.f_rhs = [](double x, double y, parameters_t p) { return p.gamma * (2.0 * x + 3.0 * y + 1.0); };
    f_lin.firstBC = { [](double x, double y) { return 2.0 * x + 3.0 * y + 1.0; } };
    f_lin.thirdBC = { 
        [](double x, double y) { return 1.0 * (3.0) + 2.0 * (2.0 * x + 3.0 * y + 1.0); }, 
        [](double x, double y) { return 1.0 * (-3.0) + 2.0 * (2.0 * x + 3.0 * y + 1.0); } 
    };

    functionsBC f_quad;
    f_quad.u_exact = [](double x, double y) { return x * x + y * y; };
    f_quad.f_rhs = [](double x, double y, parameters_t p) { return -p.lambda * 4.0 + p.gamma * (x * x + y * y); };
    f_quad.firstBC = { [](double x, double y) { return x * x + y * y; } };
    f_quad.thirdBC = { 
        [](double x, double y) { return 1.0 * (2.0 * y) + 2.0 * (x * x + y * y); }, 
        [](double x, double y) { return 1.0 * (-2.0 * y) + 2.0 * (x * x + y * y); } 
    };

    functionsBC f_trio;
    f_trio.u_exact = [](double x, double y) { return std::pow(x, 3) + std::pow(y, 3); };
    f_trio.f_rhs = [](double x, double y, parameters_t p) { 
        return -p.lambda * (6.0 * x + 6.0 * y) + p.gamma * (std::pow(x, 3) + std::pow(y, 3)); 
    };
    f_trio.firstBC = { [](double x, double y) { return std::pow(x, 3) + std::pow(y, 3); } };
    f_trio.thirdBC = { 
        [](double x, double y) { return 1.0 * (3.0 * std::pow(y, 2)) + 2.0 * (std::pow(x, 3) + std::pow(y, 3)); }, 
        [](double x, double y) { return 1.0 * (-3.0 * std::pow(y, 2)) + 2.0 * (std::pow(x, 3) + std::pow(y, 3)); } 
    };

    functionsBC f_quatro;
    f_quatro.u_exact = [](double x, double y) { return std::pow(x, 4) + std::pow(y, 4); };
    f_quatro.f_rhs = [](double x, double y, parameters_t p) { 
        return -p.lambda * (12.0 * std::pow(x, 2) + 12.0 * std::pow(y, 2)) + p.gamma * (std::pow(x, 4) + std::pow(y, 4)); 
    };
    f_quatro.firstBC = { [](double x, double y) { return std::pow(x, 4) + std::pow(y, 4); } };
    f_quatro.thirdBC = { 
        [](double x, double y) { return 1.0 * (4.0 * std::pow(y, 3)) + 2.0 * (std::pow(x, 4) + std::pow(y, 4)); }, 
        [](double x, double y) { return 1.0 * (-4.0 * std::pow(y, 3)) + 2.0 * (std::pow(x, 4) + std::pow(y, 4)); } 
    };

    functionsBC f_cinco;
    f_cinco.u_exact = [](double x, double y) { return std::pow(x, 5) + std::pow(y, 5); };
    f_cinco.f_rhs = [](double x, double y, parameters_t p) { 
        return -p.lambda * (20.0 * std::pow(x, 3) + 20.0 * std::pow(y, 3)) + p.gamma * (std::pow(x, 5) + std::pow(y, 5)); 
    };
    f_cinco.firstBC = { [](double x, double y) { return std::pow(x, 5) + std::pow(y, 5); } };
    f_cinco.thirdBC = { 
        [](double x, double y) { return 1.0 * (5.0 * std::pow(y, 4)) + 2.0 * (std::pow(x, 5) + std::pow(y, 5)); }, 
        [](double x, double y) { return 1.0 * (-5.0 * std::pow(y, 4)) + 2.0 * (std::pow(x, 5) + std::pow(y, 5)); } 
    };

    functionsBC f_seis;
    f_seis.u_exact = [](double x, double y) { return std::pow(x, 6) + std::pow(y, 6); };
    f_seis.f_rhs = [](double x, double y, parameters_t p) { 
        return -p.lambda * (30.0 * std::pow(x, 4) + 30.0 * std::pow(y, 4)) + p.gamma * (std::pow(x, 6) + std::pow(y, 6)); 
    };
    f_seis.firstBC = { [](double x, double y) { return std::pow(x, 6) + std::pow(y, 6); } };
    f_seis.thirdBC = { 
        [](double x, double y) { return 1.0 * (6.0 * std::pow(y, 5)) + 2.0 * (std::pow(x, 6) + std::pow(y, 6)); }, 
        [](double x, double y) { return 1.0 * (-6.0 * std::pow(y, 5)) + 2.0 * (std::pow(x, 6) + std::pow(y, 6)); } 
    };

    functionsBC f_sin;
    f_sin.u_exact = [](double x, double y) { return sin(x + y); };
    f_sin.f_rhs = [](double x, double y, parameters_t p) { return p.lambda * 2.0 * sin(x + y) + p.gamma * sin(x + y); };
    f_sin.firstBC = { [](double x, double y) { return sin(x + y); } };
    f_sin.thirdBC = { 
        [](double x, double y) { return 1.0 * (cos(x + y)) + 2.0 * (sin(x + y)); }, 
        [](double x, double y) { return 1.0 * (-cos(x + y)) + 2.0 * (sin(x + y)); } 
    };

    std::cout << "\n----------------------------------------------------------------------\n";
    std::cout << "ИССЛЕДОВАНИЕ ПОРЯДКА АППРОКСИМАЦИИ (Сетки N = 10, 20, 40)\n";
    std::cout << "----------------------------------------------------------------------\n";
    
    // Запускаем расчет сходимости для всех функций
    run_convergence_study("1. Linear", f_lin, area, params, eps, maxIter);
    run_convergence_study("2. Quad", f_quad, area, params, eps, maxIter);
    run_convergence_study("3. Cubic", f_trio, area, params, eps, maxIter);
    run_convergence_study("4. Power_4", f_quatro, area, params, eps, maxIter);
    run_convergence_study("5. Power_5", f_cinco, area, params, eps, maxIter);
    run_convergence_study("6. Power_6", f_seis, area, params, eps, maxIter);
    run_convergence_study("7. Sin(x+y)", f_sin, area, params, eps, maxIter);

    std::cout << "----------------------------------------------------------------------\n";
    return 0;
}