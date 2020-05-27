#include <bits/stdc++.h>
using namespace std;

#include "Polygon.h"

#define all(s) s.begin(), s.end()
#define pb push_back
#define ii pair<int, int>
#define x first
#define y second
#define bit(x, y) ((x >> y) & 1)

int main() {
    ios::sync_with_stdio(false); cin.tie(0); 
    cout.tie(0);
    freopen("in.txt", "r", stdin); 
    freopen("out.txt", "w", stdout); 
    vector<Point> circleA;
    for (int i = 0; i < 100; i++) {
        circleA.push_back(Point(1, 0).rotate(i * PI / 50));
    }

    vector<Point> circleB;
    Point bias = Point(2, 0);
    for (int i = 0; i < 100; i++) {
        circleB.push_back(bias + Point(1.5, 0).rotate(i * PI / 50));
    }

    vector<Point> intersection = ConvexPolygon::convex_intersect(circleA, circleB);
    
    ConvexPolygon cvp = ConvexPolygon(intersection);

    cvp.extend_polygon(0.3);

    cvp.printPolygon();
}
