#include "fn.h"

bool isDummyNode(area_t &area, int i, int j) {
    int L = area.elems.size();
    for (int k = 0; k < L; k++) {
        int mx0 = area.x_indices[area.elems[k].x1], mx1 = area.x_indices[area.elems[k].x2];
        int my0 = area.y_indices[area.elems[k].y1], my1 = area.y_indices[area.elems[k].y2];
        if (i >= mx0 && i <= mx1 && j >= my0 && j <= my1) {
            return false;
        }
    }
    return true;
}

void createSLAE(slae_t &slae, grid_t &grid, area_t &area, parameters_t &params, functionsBC &func) {
    int n_x = grid.x.size(), n_y = grid.y.size();
    int n = n_x * n_y;

    slae.ig = { -n_x, -1, 1, n_x };
    slae.u1.assign(n - n_x, 0.0);
    slae.u2.assign(n - 1, 0.0);
    slae.di.assign(n, 0.0);
    slae.l4.assign(n - 1, 0.0);
    slae.l5.assign(n - n_x, 0.0);
    slae.b.assign(n, 0.0);

    // Заполнение внутренних узлов 5-точечным шаблоном
    for (int j = 1; j < n_y - 1; j++) {
        for (int i = 1; i < n_x - 1; i++) {
            if (!isDummyNode(area, i, j)) {
                int k = j * n_x + i;
                double hx1 = grid.x[i] - grid.x[i - 1];
                double hx2 = grid.x[i + 1] - grid.x[i];
                double hy1 = grid.y[j] - grid.y[j - 1];
                double hy2 = grid.y[j + 1] - grid.y[j];

                slae.di[k] = params.lambda * (2.0 / (hx1 * hx2) + 2.0 / (hy1 * hy2)) + params.gamma;
                slae.u2[k - 1] = -params.lambda * 2.0 / (hx1 * (hx1 + hx2));   // левый
                slae.l4[k]     = -params.lambda * 2.0 / (hx2 * (hx1 + hx2));   // правый
                slae.u1[k - n_x]= -params.lambda * 2.0 / (hy1 * (hy1 + hy2));  // нижний
                slae.l5[k]     = -params.lambda * 2.0 / (hy2 * (hy1 + hy2));   // верхний
                
                slae.b[k] = func.f_rhs(grid.x[i], grid.y[j], params);
            }
        }
    }
}

void applyBoundaryConds(area_t &area, grid_t &grid, slae_t &slae, parameters_t &par, functionsBC &func) {
    int n_x = grid.x.size(), n_y = grid.y.size();
    double EPS = 1e-7;

    // Геометрия Т-области из area.txt
    double cap_x_min = area.x_lines[0], cap_x_max = area.x_lines[3];
    double cap_y_min = area.y_lines[1], cap_y_max = area.y_lines[2];
    double base_x_min = area.x_lines[1], base_x_max = area.x_lines[2];
    double base_y_min = area.y_lines[0], base_y_max = area.y_lines[1];

    for (int j = 0; j < n_y; j++) {
        for (int i = 0; i < n_x; i++) {
            double x = grid.x[i], y = grid.y[j];
            int k = j * n_x + i;

            if (isDummyNode(area, i, j)) continue;

            bool top_bound = (abs(y - cap_y_max) < EPS);
            bool bottom_bound = (abs(y - base_y_min) < EPS);
            bool left_cap = (abs(x - cap_x_min) < EPS && y >= cap_y_min - EPS);
            bool right_cap = (abs(x - cap_x_max) < EPS && y >= cap_y_min - EPS);
            bool left_base = (abs(x - base_x_min) < EPS && y <= base_y_max + EPS);
            bool right_base = (abs(x - base_x_max) < EPS && y <= base_y_max + EPS);
            bool inner_corners = (abs(y - cap_y_min) < EPS && (x < base_x_min - EPS || x > base_x_max + EPS));

            if (left_cap || right_cap || left_base || right_base || inner_corners) {
                // Краевое условие 1-го рода (Дирихле)
                slae.di[k] = 1.0;
                if (k - n_x >= 0) slae.u1[k - n_x] = 0;
                if (k - 1 >= 0) slae.u2[k - 1] = 0;
                if (k < slae.l4.size()) slae.l4[k] = 0;
                if (k < slae.l5.size()) slae.l5[k] = 0;
                slae.b[k] = func.firstBC[0](x, y);
            } else if (top_bound) {
                // Краевое условие 3-го рода (Робин) на верхней границе - 2-й порядок
                double hx1 = grid.x[i] - grid.x[i - 1];
                double hx2 = grid.x[i + 1] - grid.x[i];
                double hy = grid.y[j] - grid.y[j - 1];
                
                slae.di[k] = par.lambda * (2.0 / (hx1 * hx2) + 2.0 / (hy * hy)) + par.gamma + 2.0 * par.beta / hy;
                if (k - 1 >= 0) slae.u2[k - 1] = -par.lambda * 2.0 / (hx1 * (hx1 + hx2));
                if (k < slae.l4.size()) slae.l4[k] = -par.lambda * 2.0 / (hx2 * (hx1 + hx2));
                if (k - n_x >= 0) slae.u1[k - n_x] = -par.lambda * 2.0 / (hy * hy);
                
                slae.b[k] = func.f_rhs(x, y, par) + 2.0 * func.thirdBC[0](x, y) / hy;
                if (k < slae.l5.size()) slae.l5[k] = 0;
            } else if (bottom_bound) {
                // Краевое условие 3-го рода (Робин) на нижней границе - 2-й порядок
                double hx1 = grid.x[i] - grid.x[i - 1];
                double hx2 = grid.x[i + 1] - grid.x[i];
                double hy = grid.y[j + 1] - grid.y[j];
                
                slae.di[k] = par.lambda * (2.0 / (hx1 * hx2) + 2.0 / (hy * hy)) + par.gamma + 2.0 * par.beta / hy;
                if (k - 1 >= 0) slae.u2[k - 1] = -par.lambda * 2.0 / (hx1 * (hx1 + hx2));
                if (k < slae.l4.size()) slae.l4[k] = -par.lambda * 2.0 / (hx2 * (hx1 + hx2));
                if (k < slae.l5.size()) slae.l5[k] = -par.lambda * 2.0 / (hy * hy);
                
                slae.b[k] = func.f_rhs(x, y, par) + 2.0 * (func.thirdBC.size() > 1 ? func.thirdBC[1](x, y) : func.thirdBC[0](x, y)) / hy;
                if (k - n_x >= 0) slae.u1[k - n_x] = 0;
            }
        }
    }
}

void clearDummyNodes(slae_t &slae, area_t &area, int n_x, int n_y) {
    for (int j = 0; j < n_y; j++) {
        for (int i = 0; i < n_x; i++) {
            if (isDummyNode(area, i, j)) {
                int k = j * n_x + i;
                slae.di[k] = 1.0;
                if (k - n_x >= 0) slae.u1[k - n_x] = 0;
                if (k - 1 >= 0) slae.u2[k - 1] = 0;
                if (k < slae.l4.size()) slae.l4[k] = 0;
                if (k < slae.l5.size()) slae.l5[k] = 0;
                slae.b[k] = 0.0;
            }
        }
    }
}