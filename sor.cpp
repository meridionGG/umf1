#include "fn.h"
#include <iomanip>

double calcNorm(const std::vector<double> &b) {
    double sum = 0;
    for (double val : b) sum += val * val;
    return sqrt(sum);
}

double calcDiscrepancy(const slae_t &slae, const std::vector<double> &u) {
    double res = 0;
    int n = slae.di.size();
    for (int i = 0; i < n; i++) {
        double sum = slae.di[i] * u[i];
        if (i + slae.ig[0] >= 0) sum += slae.u1[i + slae.ig[0]] * u[i + slae.ig[0]];
        if (i + slae.ig[1] >= 0) sum += slae.u2[i + slae.ig[1]] * u[i + slae.ig[1]];
        if (i + slae.ig[2] < n) sum += slae.l4[i] * u[i + slae.ig[2]];
        if (i + slae.ig[3] < n) sum += slae.l5[i] * u[i + slae.ig[3]];
        res += (slae.b[i] - sum) * (slae.b[i] - sum);
    }
    return res;
}

// Алгоритм прогонки Томаса
void tridiagonal_solve(const std::vector<double>& a, const std::vector<double>& c, 
                       const std::vector<double>& b, const std::vector<double>& f, std::vector<double>& x) {
    int n = f.size();
    std::vector<double> alpha(n, 0.0), beta(n, 0.0);
    alpha[0] = -b[0] / c[0];
    beta[0] = f[0] / c[0];
    for (int i = 1; i < n - 1; ++i) {
        double denom = c[i] + a[i] * alpha[i - 1];
        alpha[i] = -b[i] / denom;
        beta[i] = (f[i] - a[i] * beta[i - 1]) / denom;
    }
    x[n - 1] = (f[n - 1] - a[n - 1] * beta[n - 2]) / (c[n - 1] + a[n - 1] * alpha[n - 2]);
    for (int i = n - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }
}

// Метод блочной релаксации (по строкам)
void method_BlockRelaxation(slae_t &slae, std::vector<double> &u, std::vector<double> &uk, double eps, int maxIter, grid_t &grid) {
    int n_x = grid.x.size(), n_y = grid.y.size();
    int n = slae.di.size();
    
    u.assign(n, 0.0);
    uk = u;
    
    double norm = calcNorm(slae.b);
    double res = 1.0;
    int iter = 0;
    
    while (res > eps && iter < maxIter) {
        for (int j = 0; j < n_y; ++j) {
            std::vector<double> a(n_x, 0.0), c(n_x, 0.0), b_sub(n_x, 0.0), f(n_x, 0.0), x_line(n_x, 0.0);
            
            for (int i = 0; i < n_x; ++i) {
                int k = j * n_x + i;
                c[i] = slae.di[k];
                if (i > 0) a[i] = slae.u2[k - 1]; 
                if (i < n_x - 1) b_sub[i] = slae.l4[k];
                
                f[i] = slae.b[k];
                // Учитываем связи с соседними строками
                if (j > 0) f[i] -= slae.u1[k - n_x] * u[k - n_x];
                if (j < n_y - 1) f[i] -= slae.l5[k] * uk[k + n_x]; 
            }
            
            tridiagonal_solve(a, c, b_sub, f, x_line);
            
            for (int i = 0; i < n_x; ++i) {
                u[j * n_x + i] = x_line[i];
            }
        }
        
        res = sqrt(calcDiscrepancy(slae, u) / norm);
        uk = u;
        iter++;
    }
    std::cout << "Сходимость достигнута за " << iter << " итераций.\nНевязка: " << std::scientific << res << "\n";
}

void output(area_t &area, std::vector<double> &u, grid_t &g, functionsBC &func, std::string filename) {
    auto &x = g.x, &y = g.y;
    const int n_x = x.size(), n_y = y.size();
    
    std::ofstream fout(filename);
    if (!fout.is_open()) {
        std::cerr << "Ошибка: Не удалось создать файл " << filename << "!\n";
        return;
    }

    fout << "N\tx\ty\tu\tu*\t|u - u*|\n";
    
    for (int j = 0; j < n_y; j++) {
        for (int i = 0; i < n_x; i++) {
            int k = j * n_x + i;
            if (!isDummyNode(area, i, j)) {
                double exact = func.u_exact(x[i], y[j]);
                double err = std::abs(u[k] - exact);
                fout << k << "\t" << x[i] << "\t" << y[j] << "\t" 
                     << std::scientific << u[k] << "\t" << exact << "\t" << err << "\n";
            }
        }
    }
    fout.close();
}