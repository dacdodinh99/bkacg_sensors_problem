#define EPS 1e-6
const double PI = acos(-1.0);

double DEG_to_RAD(double d) { return d * PI / 180.0; }
double RAD_to_DEG(double r) { return r * 180.0 / PI; }

inline int cmp(double a, double b) {
    return (a < b - EPS) ? -1 : ((a > b + EPS) ? 1 : 0);
}

struct Point {
    double x, y;
    Point() { x = y = 0.0; }
    Point(double x, double y) : x(x), y(y) {}

    Point operator + (const Point& a) const { return Point(x+a.x, y+a.y); }
    Point operator - (const Point& a) const { return Point(x-a.x, y-a.y); }
    Point operator * (double k) const { return Point(x*k, y*k); }
    Point operator / (double k) const { return Point(x/k, y/k); }

    double operator * (const Point& a) const { return x*a.x + y*a.y; } // dot product
    double operator % (const Point& a) const { return x*a.y - y*a.x; } // cross product

    int cmp(Point q) const { if (int t = ::cmp(x,q.x)) return t; return ::cmp(y,q.y); }

    #define Comp(x) bool operator x (Point q) const { return cmp(q) x 0; }
    Comp(>) Comp(<) Comp(==) Comp(>=) Comp(<=) Comp(!=)
    #undef Comp

    Point conj() { return Point(x, -y); }
    double norm() { return x*x + y*y; }

    // Note: There are 2 ways for implementing len():
    // 1. sqrt(norm()) --> fast, but inaccurate (produce some values that are of order X^2)
    // 2. hypot(x, y) --> slow, but much more accurate
    double len() { return sqrt(norm()); }

    Point rotate(double alpha) {
        double cosa = cos(alpha), sina = sin(alpha);
        return Point(x * cosa - y * sina, x * sina + y * cosa);
    }
};

struct ConvexPolygon { // The order for the vertices of all ConvexPolygon P, Q has to be counter-clockwise
    const int MAXN = 100;
    typedef vector< Point > Polygon;

public:
    ConvexPolygon(Polygon _polygon) {
        this->polygon = _polygon;
    }

    void extend_polygon(double r) {
        assert(polygon.size() > 2);

        int sz = polygon.size();
        vector<Point> temp_res;
        
        for (int i = 0; i < sz; i++) {
            Point u = polygon[i];
            Point v = polygon[(i + 1) % sz];

            // Translate uv
            Point uv = v - u;
            Point normal = Point(-uv.y, uv.x);
            Point test_dir = v + normal;
            if (ccw(u, v, test_dir)) normal = Point(0, 0) - normal;
            normal = normal * r / normal.len();
            temp_res.push_back(u + normal); 
            temp_res.push_back(v + normal);
        }

        vector<Point> final_res;
        
        sz = temp_res.size();
        for (int i = 0; i < temp_res.size(); i++) {
            final_res.push_back(temp_res[i]);

            if (i % 2 == 0) continue;

            Point origin = polygon[((i + 1) / 2) % polygon.size()];

            double angle = getAngle(origin, temp_res[i], temp_res[(i + 1) % sz]);

            Point pivot = temp_res[i] - origin;

            for (double j = PI / 50; j < angle - EPS; j += PI / 50) {
                Point delta = pivot.rotate(j);
                final_res.push_back(origin + delta);
            }

        }

        this->polygon = final_res; 

        simplifyPolygon();

    }

    Polygon getPolygon() {
        return this->polygon;
    }


    // Execute[{"A = (1, 2)"," B = (3, 4)"," C = (5, 6)"}]
    void printPolygon() {
        cout << polygon.size() << endl;
        for (int i = 0; i < polygon.size(); i++) {
            cout << polygon[i].x << " " << polygon[i].y;
            if (i + 1 != polygon.size()) cout << endl; 
        }
    }

    static Polygon convex_intersect(Polygon P, Polygon Q) {
        const int n = P.size(), m = Q.size();
        int a = 0, b = 0, aa = 0, ba = 0;
        enum { Pin, Qin, Unknown } in = Unknown;
        Polygon R;
        do {
            int a1 = (a+n-1) % n, b1 = (b+m-1) % m;
            double C = (P[a] - P[a1]) % (Q[b] - Q[b1]);
            double A = (P[a1] - Q[b]) % (P[a] - Q[b]);
            double B = (Q[b1] - P[a]) % (Q[b] - P[a]);
            Point r;
            if (intersect_1pt(P[a1], P[a], Q[b1], Q[b], r)) {
                if (in == Unknown) aa = ba = 0;
                R.push_back( r );
                in = B > 0 ? Pin : A > 0 ? Qin : in;
            }
            if (C == 0 && B == 0 && A == 0) {
                if (in == Pin) { b = (b + 1) % m; ++ba; }
                else            { a = (a + 1) % m; ++aa; }
            } else if (C >= 0) {
                if (A > 0) { if (in == Pin) R.push_back(P[a]); a = (a+1)%n; ++aa; }
                else       { if (in == Qin) R.push_back(Q[b]); b = (b+1)%m; ++ba; }
            } else {
                if (B > 0) { if (in == Qin) R.push_back(Q[b]); b = (b+1)%m; ++ba; }
                else       { if (in == Pin) R.push_back(P[a]); a = (a+1)%n; ++aa; }
            }
        } while ( (aa < n || ba < m) && aa < 2*n && ba < 2*m );
        if (in == Unknown) {
            if (in_convex(Q, P[0])) return P;
            if (in_convex(P, Q[0])) return Q;
        }
        return R;
    }

private:
    Polygon polygon;

    static int ccw(Point a, Point b, Point c) {
        return cmp((b-a)%(c-a),0);
    }

    static double getAngle(Point a, Point b, Point c) {
        Point ab = b - a;
        Point ac = c - a;
        return acos((ab * ac / (ab.len() * ac.len())));
    }

    static bool intersect_1pt(Point a, Point b,
        Point c, Point d, Point &r) {
        double D =  (b - a) % (d - c);
        if (cmp(D, 0) == 0) return false;
        double t =  ((c - a) % (d - c)) / D;
        double s = -((a - c) % (b - a)) / D;
        r = a + (b - a) * t;
        return cmp(t, 0) >= 0 && cmp(t, 1) <= 0 && cmp(s, 0) >= 0 && cmp(s, 1) <= 0;
    }

    #define Det(a,b,c) ((double)(b.x-a.x)*(double)(c.y-a.y)-(double)(b.y-a.y)*(c.x-a.x))
    static bool in_convex(vector<Point>& l, Point p){
        int a = 1, b = l.size()-1, c;
        if (Det(l[0], l[a], l[b]) > 0) swap(a,b);
        // Allow on edge --> if (Det... > 0 || Det ... < 0)
        if (Det(l[0], l[a], p) >= 0 || Det(l[0], l[b], p) <= 0) return false;
        while(abs(a-b) > 1) {
            c = (a+b)/2;
            if (Det(l[0], l[c], p) > 0) b = c; else a = c;
        }
        // Alow on edge --> return Det... <= 0
        return Det(l[a], l[b], p) < 0;
    }

    double getArea(Point a, Point b, Point c) {
        return abs((b - a) % (c - a));
    }

    void simplifyPolygon() {
        while (polygon.size() > MAXN) {
            double minArea = 1e18;
            int id = 0;
            int sz = polygon.size();

            for (int i = 0; i < sz; i++) {
                Point u = polygon[i];
                Point v = polygon[(i + sz - 1) % sz];
                Point t = polygon[(i + 1) % sz];

                double area = getArea(u, v, t);
                if (area < minArea) {
                    minArea = area;
                    id = i;
                }
            }

            polygon.erase(polygon.begin() + id);
        }
    }

};




