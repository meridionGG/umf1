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

// Добавили std::string out_filename
double run_test(functionsBC func, area_t &area, parameters_t &params, int N, double eps, int maxIter, std::string out_filename) {
    write_grid(N);
    grid_t grid;
    read_grid("in/grid.txt", area, grid);

    slae_t slae;
    std::vector<double> u, uk;

    createSLAE(slae, grid, area, params, func);
    applyBoundaryConds(area, grid, slae, params, func);
    clearDummyNodes(slae, area, grid.x.size(), grid.y.size());

    // Решаем СЛАУ
    method_BlockRelaxation(slae, u, uk, eps, maxIter, grid);

    // СОХРАНЯЕМ ДАННЫЕ В ФАЙЛ
    output(area, u, grid, func, out_filename);

    // Считаем максимальную ошибку для вывода в консоль
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

int main() {
    const double eps = 1e-13;
    const int maxIter = 50000;

    area_t area;
    parameters_t params;
    
    std::cout << "Чтение базовых параметров !!!!!!..." << std::endl;
    read_area("in/area.txt", area);
    read_params("in/params.txt", params);
    
    // --- ЖЕЛЕЗОБЕТОННЫЙ КОСТЫЛЬ ---
    // Принудительно задаем параметры (на случай, если params.txt читается криво)
    params.lambda = 1.0;
    params.gamma = 1.0;
    params.beta = 2.0;

    // Проверяем, прочиталась ли вообще геометрия
    if (area.x_lines.size() < 4 || area.y_lines.size() < 3) {
        std::cerr << "КРИТИЧЕСКАЯ ОШИБКА: Файл in/area.txt не прочитался или пуст!\n";
        std::cerr << "Убедитесь, что папка 'in' лежит ровно там же, откуда запускается программа.\n";
        return 1; // Экстренно завершаем программу, чтобы не получить -nan
    }
    // -------------------------------
    
    std::cout << "Запуск автоматических тестов...\n" << std::endl;

    // --- Линейный полином ---
    functionsBC f_lin;
    f_lin.u_exact = [](double x, double y) { return 2.0 * x + 3.0 * y + 1.0; };
    f_lin.f_rhs = [](double x, double y, parameters_t p) { return p.gamma * (2.0 * x + 3.0 * y + 1.0); };
    f_lin.firstBC = { [](double x, double y) { return 2.0 * x + 3.0 * y + 1.0; } };
    f_lin.thirdBC = { [](double x, double y) { return 2.0 * x + 3.0 * y + 1.0; } };

    // --- Квадратичный полином ---
    functionsBC f_quad;
    f_quad.u_exact = [](double x, double y) { return x * x + y * y; };
    f_quad.f_rhs = [](double x, double y, parameters_t p) { return -p.lambda * 4.0 + p.gamma * (x * x + y * y); };
    f_quad.firstBC = { [](double x, double y) { return x * x + y * y; } };
    f_quad.thirdBC = { [](double x, double y) { return x * x + y * y; } };

    // --- Синус ---
    functionsBC f_sin;
    f_sin.u_exact = [](double x, double y) { return sin(x + y); };
    f_sin.f_rhs = [](double x, double y, parameters_t p) { return p.lambda * 2.0 * sin(x + y) + p.gamma * sin(x + y); };
    f_sin.firstBC = { [](double x, double y) { return sin(x + y); } };
    f_sin.thirdBC = { [](double x, double y) { return sin(x + y); } };


    std::cout << "---------------------------------------------------------\n";
    std::cout << "ТАБЛИЦА 1: Тестирование на полиномах\n";
    std::cout << "---------------------------------------------------------\n";
    
    // ПЕРЕДАЕМ ИМЕНА ФАЙЛОВ
    double err_lin = run_test(f_lin, area, params, 10, eps, maxIter, "out/out_linear.txt");
    std::cout << std::left << std::setw(30) << "1. Линейный:" << std::scientific << err_lin << "\n";

    double err_quad = run_test(f_quad, area, params, 10, eps, maxIter, "out/out_quad.txt");
    std::cout << std::left << std::setw(30) << "2. Квадратичный:" << std::scientific << err_quad << "\n";


    std::cout << "\n---------------------------------------------------------\n";
    std::cout << "ТАБЛИЦА 2: Исследование порядка аппроксимации ( u = sin(x+y) )\n";
    std::cout << "---------------------------------------------------------\n";
    
    std::vector<int> Ns = {10, 20, 40};
    std::vector<double> errs(3);

    for (size_t i = 0; i < Ns.size(); ++i) {
        // Динамическое имя файла для каждой сетки синуса
        std::string filename = "out/out_sin_" + std::to_string(Ns[i]) + ".txt";
        errs[i] = run_test(f_sin, area, params, Ns[i], eps, maxIter, filename);
        
        std::cout << "Узлов (N) = " << std::setw(2) << Ns[i] << " | Ошибка Eh = " << std::fixed << std::setprecision(6) << errs[i];
        
        if (i > 0) std::cout << " | Порядок p = " << std::setprecision(3) << log2(errs[i - 1] / errs[i]) << "\n";
        else std::cout << " | Порядок p = -\n";
    }
    
    std::cout << "---------------------------------------------------------\n";
    std::cout << "Готово! Проверьте папку 'out' - там появилось 5 файлов с результатами." << std::endl;
    return 0;
}