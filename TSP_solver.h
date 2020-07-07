#pragma once


#include <bits/stdc++.h>
using namespace std;
#include "Polygon.h"

#define all(s) s.begin(), s.end()
#define pb push_back
#define ii pair<int, int>
#define bit(x, y) ((x >> y) & 1)

const int N = 15;

vector<double> dp[1 << N][N];
vector<Point> polygons[N];

double TSP_solver(const vector<ConvexPolygon> &a) {
    int n = (int) a.size();
    if (n <= 1) return 0.0;

    for (int i = 0; i < n; i++) {
        polygons[i] = a[i].getPolygon();
        // cout << "Da giac thu:" << i << endl;
        // for (auto e : polygons[i]) cout << e.x << ' ' << e.y << endl;
    }

    int id = 0;
    for (int i = 1; i < n; i++) {
        if (polygons[i].size() < polygons[id].size()) id = i;
    }
    swap(polygons[id], polygons[0]);

    double res = 1e18;

    for (int start = 0; start < 1; start++) {
        for (int mask = 0; mask < (1 << n); mask++) {
            for (int i = 0; i < n; i++) dp[mask][i].resize((int) polygons[i].size(), 1e18);
        }

        dp[1][0][start] = 0;

        for (int mask = 2; mask < (1 << n); mask++) {
            for (int i = 1; i < n; i++) {
                if (bit(mask, i) == 0) continue;
                for (int u = 0; u < (int) polygons[i].size(); u++) {
                    for (int j = 0; j < n; j++) {
                        if (j == i || bit(mask, j) == 0) continue;
                        for (int v = 0; v < (int) polygons[j].size(); v++) {
                            double len = (polygons[i][u] - polygons[j][v]).len();
                            double pre = dp[mask][i][u];
                            dp[mask][i][u] = min(dp[mask][i][u], dp[mask ^ (1 << i)][j][v] + len);
                            // if (dp[mask][i][u] < 1e-6) {
                            //     cout << polygons[i][u].x << ' ' << polygons[i][u].y << endl;
                            //     cout << polygons[j][v].x << ' ' << polygons[j][v].y << endl;
                            // }
                        }
                    }
                }
            }
        }

        double tmp_res = 1e18;
        for (int i = 1; i < n; i++) {
            for (int j = 0; j < (int) polygons[i].size(); j++) {
                double len = (polygons[i][j] - polygons[0][start]).len();
                tmp_res = min(tmp_res, len + dp[(1 << n) - 1][i][j]);
            }
        }

        res = min(res, tmp_res);
    }
    // cout << "Dap so TSP: " << res << endl;
    return res;
}

// int main() {
//     ios::sync_with_stdio(false); cin.tie(0); 
//     cout.tie(0);
//     freopen("in.txt", "r", stdin); 
//     int n;
//     cin >> n;
//     vector<ConvexPolygon> polygons;
//     for (int i = 1; i <= n; i++) {
//         int nVer; cin >> nVer;
//         vector<Point> polygon;
//         while (nVer--) {
//             int x, y;
//             cin >> x >> y;
//             polygon.pb(Point(x, y));
//         }
//         polygons.pb(ConvexPolygon(polygon));
//     }

//     cout << TSP_solver(polygons);
//     return 0;
// }