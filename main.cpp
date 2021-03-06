#include <bits/stdc++.h>
using namespace std;
#include "Polygon.h"
#include "TSP_solver.h"


#define all(s) s.begin(), s.end()
#define pb push_back
#define bit(x, y) ((x >> y) & 1)

double extend(int k, vector<ConvexPolygon> polygons, double R) {
    
    polygons[k].extend_polygon(R);

    return TSP_solver(polygons); 
}

double hop(int x, int y, vector<ConvexPolygon> polygons, double r) {
    if (x > y) swap(x, y);
    vector<ConvexPolygon> new_polygons;
    for (int i = 0; i < (int) polygons.size(); i++) {
        if (i != x && i != y) new_polygons.push_back(polygons[i]);
    }

    ConvexPolygon new_polygon = ConvexPolygon::convex_intersect(polygons[x].getPolygon(), polygons[y].getPolygon());

    if (new_polygon.getPolygon().size() == 0) return 1e18;

    new_polygon.extend_polygon(r);

    new_polygons.push_back(new_polygon);

    return TSP_solver(new_polygons);
}

int main() {
    ios::sync_with_stdio(false); cin.tie(0); 
    cout.tie(0);
    freopen("in.txt", "r", stdin); 
    int N;
    double R, L;

    cin >> N >> R >> L;
    vector<ConvexPolygon> polygons;
    for (int i = 1; i <= N; i++) {
        double x, y;
        cin >> x >> y;
        ConvexPolygon new_polygon = ConvexPolygon(x, y, R);
        polygons.push_back(new_polygon);
    }

    int ans = 0;

    while (TSP_solver(polygons) > L) {
        cout << "So sensor hien tai: " << ans << endl;
        pair<int, int> mem = {-1, -1};
        double best = 1e9;
        
        // Check extend polygons
        for (int i = 0; i < (int) polygons.size(); i++) {
            double foo = extend(i, polygons, R);
            if (best > foo) {
                mem = {i, -1};
                best = foo;
            }
        }

        // Hop
        for (int i = 0; i < (int) polygons.size(); i++) {
            for (int j = i + 1; j < (int) polygons.size(); j++) {
                double foo = hop(i, j, polygons, R);
                if (best > foo) {
                    mem = {i, j};
                    best = foo;
                }
            }
        }


        if (mem.second == -1) {
            cout << "No da giac:" << mem.first << endl;
            polygons[mem.first].extend_polygon(R);
        } else {
            cout << "Hop da giac" << ' ' << mem.first << ' ' << mem.second << endl;

            vector<ConvexPolygon> new_polygons;
            for (int i = 0; i < (int) polygons.size(); i++) {
                if (i != mem.first && i != mem.second) new_polygons.push_back(polygons[i]);
            }


            ConvexPolygon new_polygon = ConvexPolygon::convex_intersect(polygons[mem.first].getPolygon(), polygons[mem.second].getPolygon());

            new_polygon.extend_polygon(R);


            new_polygons.push_back(new_polygon);

            polygons.clear();
            for (auto e : new_polygons) polygons.push_back(e);
            new_polygons.clear();

        }

        ans++;
    }
    cout << TSP_solver(polygons) << endl;
    cout << "Ket qua cuoi cung: " << ans << endl;
    return 0;
}