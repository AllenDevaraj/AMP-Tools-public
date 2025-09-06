#include "MyBugAlgorithm.h"
#include <algorithm>
#include <limits>

// Constants for numerical stability
static constexpr double EPS = 1e-9;
static constexpr double LEAVE_EPS = 1e-3; // How far to nudge away from the wall
static constexpr double IMPROVE_MIN = 1e-5; // Required improvement to leave

amp::Path2D MyBugAlgorithm::plan(const amp::Problem2D& problem) {
    const Eigen::Vector2d qs = problem.q_init;
    const Eigen::Vector2d qg = problem.q_goal;
    const auto& obstacles = problem.obstacles;

    std::vector<Eigen::Vector2d> polyline;
    polyline.push_back(qs);

    Eigen::Vector2d current_pos = qs;
    int loop_guard = 0;

    while ((current_pos - qg).norm() > 1e-6 && loop_guard++ < 200) {
        // 1. Try to move straight to goal, find the first hit
        Hit hit = findFirstHit(current_pos, qg, obstacles);
        if (!hit.hit) {
            appendSegment(polyline, current_pos, qg);
            break; // Path to goal is clear
        }
        appendSegment(polyline, current_pos, hit.point);
        current_pos = hit.point;

        // 2. Survey the entire obstacle boundary to find the optimal leave point (L*)
        Eigen::Vector2d Lstar; 
        size_t L_edgeIdx;
        surveyBoundaryForLstar(obstacles[hit.obsIdx], turn_, qg, Lstar, L_edgeIdx);

        // 3. If L* isn't meaningfully better than the hit point, we are likely trapped.
        if (dgoal(hit.point, qg) - dgoal(Lstar, qg) < IMPROVE_MIN) {
            break;
        }

        // 4. Generate the path along the boundary from the hit point to L*
        travelToLeavePoint(obstacles[hit.obsIdx], hit.edgeIdx, hit.point, turn_, Lstar, L_edgeIdx, polyline);
        current_pos = Lstar;

        // 5. Leave the obstacle by nudging away from the wall
        const auto V = vertsOf(obstacles[hit.obsIdx]);
        const Eigen::Vector2d outward_normal = getOutwardNormal(V, L_edgeIdx, Lstar);
        current_pos += outward_normal * LEAVE_EPS;
        appendSegment(polyline, Lstar, current_pos);
    }

    return amp::Path2D{polyline};
}

// ===========================================================================
// HELPER FUNCTION IMPLEMENTATIONS
// ===========================================================================

std::vector<Eigen::Vector2d> MyBugAlgorithm::vertsOf(const amp::Obstacle2D& obs) {
    return obs.verticesCW();
}

bool MyBugAlgorithm::segSegIntersect(const Eigen::Vector2d& p, const Eigen::Vector2d& q,
                                     const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                                     double& t, double& u, Eigen::Vector2d& x) {
    const Eigen::Vector2d r = q - p, s = b - a;
    const double denom = r.x()*s.y() - r.y()*s.x();
    const double num_t = (a.x()-p.x())*s.y() - (a.y()-p.y())*s.x();
    const double num_u = (a.x()-p.x())*r.y() - (a.y()-p.y())*r.x();
    if (std::abs(denom) < EPS) return false;
    t = num_t / denom; u = num_u / denom;
    if (t < -EPS || t > 1.0 + EPS || u < -EPS || u > 1.0 + EPS) return false;
    x = p + t * r;
    return true;
}

Eigen::Vector2d MyBugAlgorithm::closestPointOnSegment(const Eigen::Vector2d& A, const Eigen::Vector2d& B, const Eigen::Vector2d& P) {
    const Eigen::Vector2d AB = B - A;
    if (AB.squaredNorm() <= EPS) return A;
    double t = (P - A).dot(AB) / AB.squaredNorm();
    t = std::max(0.0, std::min(1.0, t));
    return A + t * AB;
}

MyBugAlgorithm::Hit MyBugAlgorithm::findFirstHit(const Eigen::Vector2d& s, const Eigen::Vector2d& g,
                                                 const std::vector<amp::Obstacle2D>& obstacles) {
    Hit best_hit; 
    double best_t = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < obstacles.size(); ++i) {
        auto V = vertsOf(obstacles[i]);
        for (size_t e = 0; e < V.size(); ++e) {
            double t, u; Eigen::Vector2d intersection_point;
            if (segSegIntersect(s, g, V[e], V[(e + 1) % V.size()], t, u, intersection_point)) {
                if (t >= EPS && t < best_t) {
                    best_t = t;
                    best_hit = {true, i, e, intersection_point};
                }
            }
        }
    }
    return best_hit;
}

void MyBugAlgorithm::surveyBoundaryForLstar(const amp::Obstacle2D& O, Turn turn,
                                            const Eigen::Vector2d& qg,
                                            Eigen::Vector2d& Lstar, size_t& L_edgeIdx) {
    auto V = vertsOf(O);
    if (V.empty()) return;

    Lstar = V[0];
    L_edgeIdx = 0;
    double min_dist = dgoal(Lstar, qg);

    for (size_t e = 0; e < V.size(); ++e) {
        Eigen::Vector2d C = closestPointOnSegment(V[e], V[(e + 1) % V.size()], qg);
        double d = dgoal(C, qg);
        if (d < min_dist) {
            min_dist = d;
            Lstar = C;
            L_edgeIdx = e;
        }
    }
}

void MyBugAlgorithm::travelToLeavePoint(const amp::Obstacle2D& O, size_t hit_edgeIdx,
                                        const Eigen::Vector2d& hit_point, Turn turn,
                                        const Eigen::Vector2d& Lstar, size_t L_edgeIdx,
                                        std::vector<Eigen::Vector2d>& polyline) {
    auto V = vertsOf(O);
    if (V.empty()) return;
    
    appendSegment(polyline, polyline.back(), hit_point);

    size_t current_edge = hit_edgeIdx;
    while (current_edge != L_edgeIdx) {
        size_t next_vertex_idx = (turn == RIGHT) ? (current_edge + 1) % V.size() : current_edge;
        appendSegment(polyline, polyline.back(), V[next_vertex_idx]);
        current_edge = (turn == RIGHT) ? (current_edge + 1) % V.size() : (current_edge + V.size() - 1) % V.size();
    }
    
    appendSegment(polyline, polyline.back(), Lstar);
}

Eigen::Vector2d MyBugAlgorithm::getOutwardNormal(const std::vector<Eigen::Vector2d>& V,
                                                 size_t edgeIdx, const Eigen::Vector2d& point) {
    const Eigen::Vector2d A = V[edgeIdx];
    const Eigen::Vector2d B = V[(edgeIdx + 1) % V.size()];
    Eigen::Vector2d edge_vec = (B - A).normalized();
    // For CW vertices, the outward normal is a right turn from the edge direction
    return Eigen::Vector2d(-edge_vec.y(), edge_vec.x());
}

void MyBugAlgorithm::appendSegment(std::vector<Eigen::Vector2d>& out,
                                   const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
    if (out.empty() || (out.back() - a).norm() > EPS) out.push_back(a);
    if ((out.back() - b).norm() > EPS) out.push_back(b);
}
