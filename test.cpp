#include <bits/stdc++.h>
using namespace std;
#include "Polygon.h"

#define all(s) s.begin(), s.end()
#define pb push_back
#define ii pair<int, int>
#define bit(x, y) ((x >> y) & 1)

int main() {
    ios::sync_with_stdio(false); cin.tie(0); 
    cout.tie(0);
    freopen("in.txt", "r", stdin); 
    int n, m;
    cin >> n;
    vector<Point> po1;
    for (int i = 0; i < n; i++) {
        Point p;
        cin >> p.x >> p.y;
        po1.pb(p);
    }
    cin >> n;
    vector<Point> po2;
    for (int i = 0; i < m; i++) {
        Point p;
        cin >> p.x >> p.y;
        po2.pb(p);
    }
    ConvexPolygon p3 = ConvexPolygon::convex_intersect(po1, po2);
    for (Point e : p3.getPolygon()) {
        cout << e.x << ' ' << e.y << endl;
    }
    return 0;
}