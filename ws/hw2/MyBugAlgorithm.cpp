#include "MyBugAlgorithm.h"
#include <algorithm>
#include <limits>

namespace {
    constexpr double EPS            = 1e-9;
    constexpr double IMPROVE_MIN    = 1e-4;  // require this much improvement at L* vs hit
    constexpr double LEAVE_EPS      = 2e-3;  // outward nudge size
    constexpr double GOAL_STEP_EPS  = 2e-3;  // extra push toward goal at leave
}

// Always return vertices in **CW** order (normalize CCW if needed)
std::vector<Eigen::Vector2d> MyBugAlgorithm::vertsCW(const amp::Obstacle2D& obs) {
    std::vector<Eigen::Vector2d> V = obs.verticesCW();
    if (V.size() < 3) return V;
    if (signedArea(V) > 0) std::reverse(V.begin(), V.end()); // make CW
    return V;
}

double MyBugAlgorithm::signedArea(const std::vector<Eigen::Vector2d>& V) {
    double A = 0.0; size_t n = V.size();
    for (size_t i = 0; i < n; ++i) {
        const auto& p = V[i];
        const auto& q = V[(i+1)%n];
        A += p.x()*q.y() - q.x()*p.y();
    }
    return 0.5 * A;
}

bool MyBugAlgorithm::segSegIntersect(const Eigen::Vector2d& p, const Eigen::Vector2d& q,
                                     const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                                     double& t, double& u, Eigen::Vector2d& x) {
    const Eigen::Vector2d r = q - p, s = b - a;
    const double denom = r.x()*s.y() - r.y()*s.x();
    const double num_t = (a.x()-p.x())*s.y() - (a.y()-p.y())*s.x();
    const double num_u = (a.x()-p.x())*r.y() - (a.y()-p.y())*r.x();
    if (std::abs(denom) < EPS) return false;     // parallel/collinear
    t = num_t / denom; u = num_u / denom;
    if (t < -EPS || t > 1.0 + EPS || u < -EPS || u > 1.0 + EPS) return false;
    x = p + t * r;
    return true;
}

Eigen::Vector2d MyBugAlgorithm::closestPointOnSegment(const Eigen::Vector2d& A,
                                                      const Eigen::Vector2d& B,
                                                      const Eigen::Vector2d& P,
                                                      double& t) {
    const Eigen::Vector2d AB = B - A;
    const double den = AB.squaredNorm();
    if (den <= EPS) { t = 0.0; return A; }
    double tau = (P - A).dot(AB) / den;
    tau = std::clamp(tau, 0.0, 1.0);
    t = tau;
    return A + tau * AB;
}

// For CW polygon, outward (exterior) normal of edge A->B is (-dy, dx)
static inline Eigen::Vector2d unitOrZero(const Eigen::Vector2d& v) {
    const double n = v.norm(); return (n>0)? v/n : Eigen::Vector2d(0,0);
}
Eigen::Vector2d MyBugAlgorithm::outwardNormal_CWEdge(const Eigen::Vector2d& A,
                                                     const Eigen::Vector2d& B) {
    Eigen::Vector2d d = unitOrZero(B-A);
    return Eigen::Vector2d(-d.y(), d.x());
}
Eigen::Vector2d MyBugAlgorithm::outwardNormalAtPointCW(const std::vector<Eigen::Vector2d>& V,
                                                       size_t e,
                                                       const Eigen::Vector2d& P) {
    const size_t n = V.size();
    const size_t eN = (e+1)%n;
    const Eigen::Vector2d& A = V[e];
    const Eigen::Vector2d& B = V[eN];
    const bool nearA = dist(P,A) <= 1e-7;
    const bool nearB = dist(P,B) <= 1e-7;
    if (nearA) {
        const size_t ePrev = (e+n-1)%n;
        return unitOrZero(outwardNormal_CWEdge(V[ePrev], V[(ePrev+1)%n]) +
                          outwardNormal_CWEdge(V[e],    V[(e+1)%n]));
    }
    if (nearB) {
        const size_t eNext = (e+1)%n;
        return unitOrZero(outwardNormal_CWEdge(V[e],    V[(e+1)%n]) +
                          outwardNormal_CWEdge(V[eNext],V[(eNext+1)%n]));
    }
    return outwardNormal_CWEdge(A,B);
}

// ---- first hit using METRIC min distance filter -----------------------------
MyBugAlgorithm::Hit
MyBugAlgorithm::firstHit(const Eigen::Vector2d& s, const Eigen::Vector2d& g,
                         const std::vector<amp::Obstacle2D>& obstacles,
                         double min_hit_dist) {
    Hit best; double bestT = std::numeric_limits<double>::infinity();
    const double rayLen = (g - s).norm();
    for (size_t i=0;i<obstacles.size();++i) {
        auto V = vertsCW(obstacles[i]); const size_t n = V.size();
        if (n<2) continue;
        for (size_t e=0;e<n;++e) {
            const auto& a = V[e]; const auto& b = V[(e+1)%n];
            double t,u; Eigen::Vector2d X;
            if (!segSegIntersect(s,g,a,b,t,u,X)) continue;
            if (t < 0.0 || t > 1.0) continue;
            const double hitDist = t * rayLen;
            if (hitDist <= min_hit_dist) continue;  // <---- new robust filter
            if (t < bestT) { bestT=t; best={true,i,e,t,X}; }
        }
    }
    return best;
}

// ---- full survey: find L* ---------------------------------------------------
void MyBugAlgorithm::surveyBoundary_Bug1(const amp::Obstacle2D& O, size_t startEdge,
                                         const Eigen::Vector2d& hitPoint, Turn turn,
                                         const Eigen::Vector2d& goal,
                                         Eigen::Vector2d& Lstar, size_t& L_edge) {
    auto V = vertsCW(O); const size_t n = V.size();
    if (n==0) { Lstar=hitPoint; L_edge=startEdge; return; }

    // With CW verts, CCW walking (LEFT) means decreasing edge index.
    auto nextEdge = [&](size_t e)->size_t {
        return (turn==LEFT) ? ((e+n-1)%n) : ((e+1)%n);
    };

    Lstar = hitPoint; L_edge = startEdge;
    double best = dist(Lstar, goal);

    // One full loop; check continuous closest point on each edge.
    size_t e = startEdge;
    for (size_t k=0;k<n;++k) {
        const auto& A = V[e]; const auto& B = V[(e+1)%n];
        double tau=0.0; Eigen::Vector2d C = closestPointOnSegment(A,B,goal,tau);
        const double d = dist(C,goal);
        if (d + 1e-12 < best - 1e-15) { best=d; Lstar=C; L_edge=e; }
        e = nextEdge(e);
    }
}

// ---- walk boundary from hit to L* -------------------------------------------
void MyBugAlgorithm::walkToLstar(const amp::Obstacle2D& O, size_t startEdge,
                                 const Eigen::Vector2d& hitPoint, Turn turn,
                                 const Eigen::Vector2d& Lstar, size_t L_edge,
                                 std::vector<Eigen::Vector2d>& polyline) {
    auto V = vertsCW(O); const size_t n = V.size(); if (n==0) return;

    auto nextEdge = [&](size_t e)->size_t {
        return (turn==LEFT) ? ((e+n-1)%n) : ((e+1)%n);
    };
    auto endVertex = [&](size_t e)->size_t {
        // end vertex for chosen traversal direction
        return (turn==LEFT) ? e : ((e+1)%n);
    };

    if (polyline.empty() || dist(polyline.back(), hitPoint) > 1e-7)
        polyline.push_back(hitPoint);

    // finish current edge to its end vertex
    size_t e = startEdge;
    size_t vEnd = endVertex(e);
    if (dist(polyline.back(), V[vEnd]) > 1e-9) polyline.push_back(V[vEnd]);

    // walk edge-by-edge until we reach L_edge
    e = nextEdge(e);
    int guard=0;
    while (e != L_edge && guard++ < static_cast<int>(2*n)+10) {
        size_t v = endVertex(e);
        if (dist(polyline.back(), V[v]) > 1e-9) polyline.push_back(V[v]);
        e = nextEdge(e);
    }

    // append L*
    if (dist(polyline.back(), Lstar) > 1e-9) polyline.push_back(Lstar);
}

// ---- append segment ----------------------------------------------------------
void MyBugAlgorithm::appendSeg(std::vector<Eigen::Vector2d>& out,
                               const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
    if (out.empty() || dist(out.back(), a) > 1e-7) out.push_back(a);
    if (dist(out.back(), b) > 1e-9) out.push_back(b);
}

// === main planner ============================================================
amp::Path2D MyBugAlgorithm::plan(const amp::Problem2D& prob) {
    const Eigen::Vector2d qs = prob.q_init;
    const Eigen::Vector2d qg = prob.q_goal;
    const auto& obstacles     = prob.obstacles;

    std::vector<Eigen::Vector2d> poly;
    poly.reserve(512);
    poly.push_back(qs);

    Eigen::Vector2d cur = qs;
    double last_d = dist(cur, qg);
    int outer_guard = 0, stagnant = 0;

    while (dist(cur, qg) > 1e-6 && outer_guard++ < 5000) {
        // Metric min-hit-dist tied to the leave nudge:
        const double MIN_HIT_DIST = 5.0 * LEAVE_EPS;

        // 1) straight cast
        Hit H = firstHit(cur, qg, obstacles, MIN_HIT_DIST);
        if (!H.hit) {
            appendSeg(poly, cur, qg);
            cur = qg;
            break;
        }
        appendSeg(poly, cur, H.point);

        // 2) full circumnavigation => L*
        Eigen::Vector2d Lstar; size_t L_edge=H.edgeIdx;
        surveyBoundary_Bug1(obstacles[H.obsIdx], H.edgeIdx, H.point, turn_, qg, Lstar, L_edge);

        // Require meaningful improvement vs hit point
        if (dist(H.point, qg) - dist(Lstar, qg) < IMPROVE_MIN) {
            // no progress possible around this obstacle
            break;
        }

        // 3) walk boundary back to L*
        walkToLstar(obstacles[H.obsIdx], H.edgeIdx, H.point, turn_, Lstar, L_edge, poly);

        // 4) peel off: outward + slight goal bias, then re-cast next loop
        const auto Vcw    = vertsCW(obstacles[H.obsIdx]);
        const auto n_out  = outwardNormalAtPointCW(Vcw, L_edge, Lstar);
        const auto g_dir  = unitOrZero(qg - Lstar);
        const Eigen::Vector2d leavePt = Lstar + LEAVE_EPS*n_out + GOAL_STEP_EPS*g_dir;

        appendSeg(poly, Lstar, leavePt);
        cur = leavePt;

        // progress guard (avoid livelock)
        const double dnow = dist(cur, qg);
        if (last_d - dnow < 1e-4) {
            if (++stagnant > 50) break;
        } else {
            stagnant = 0;
        }
        last_d = dnow;
    }

    amp::Path2D path;
    path.waypoints = std::move(poly);
    path.valid = (dist(cur, qg) <= 1e-6);
    return path;
}
