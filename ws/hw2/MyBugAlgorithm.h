#pragma once
#include <vector>
#include <Eigen/Dense>
#include <hw/HW2.h>
#include <tools/Environment.h> // amp::Problem2D, amp::Obstacle2D

// BUG-1 (left-turning / "right-hand on wall" viewpoint: keep obstacle on LEFT)
class MyBugAlgorithm : public amp::BugAlgorithm {
public:
    enum Algo { BUG1 };
    enum Turn { LEFT, RIGHT };

    explicit MyBugAlgorithm(Algo algo = BUG1, Turn turn = LEFT)
        : algo_(algo), turn_(turn) {}

    // === AMP interface ===
    amp::Path2D plan(const amp::Problem2D& prob) override;
    const char* name() const { return "MyBugAlgorithm (BUG-1 LEFT)"; }

private:
    // ---- geometry helpers ----
    static std::vector<Eigen::Vector2d> vertsCW(const amp::Obstacle2D& obs);  // always CW
    static bool segSegIntersect(const Eigen::Vector2d& p, const Eigen::Vector2d& q,
                                const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                                double& t, double& u, Eigen::Vector2d& x);
    static inline double dist(const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
        return (a - b).norm();
    }
    static double signedArea(const std::vector<Eigen::Vector2d>& V);
    static Eigen::Vector2d closestPointOnSegment(const Eigen::Vector2d& A,
                                                 const Eigen::Vector2d& B,
                                                 const Eigen::Vector2d& P,
                                                 double& t);
    static Eigen::Vector2d outwardNormal_CWEdge(const Eigen::Vector2d& A,
                                                const Eigen::Vector2d& B);
    static Eigen::Vector2d outwardNormalAtPointCW(const std::vector<Eigen::Vector2d>& V,
                                                  size_t edgeIdx,
                                                  const Eigen::Vector2d& P);

    struct Hit {
        bool hit = false;
        size_t obsIdx = 0;     // which obstacle
        size_t edgeIdx = 0;    // which edge on that obstacle (edge = [V[e], V[e+1]])
        double tAlong = 0.0;   // param along s->g
        Eigen::Vector2d point; // intersection point
    };

    // Find first hit from s→g; ignore intersections within min_hit_dist of s
    static Hit firstHit(const Eigen::Vector2d& s, const Eigen::Vector2d& g,
                        const std::vector<amp::Obstacle2D>& obstacles,
                        double min_hit_dist);

    // Full circumnavigation to find L* (global closest-to-goal boundary point)
    static void surveyBoundary_Bug1(const amp::Obstacle2D& O, size_t startEdge,
                                    const Eigen::Vector2d& hitPoint, Turn turn,
                                    const Eigen::Vector2d& goal,
                                    Eigen::Vector2d& Lstar, size_t& L_edge);

    // Walk boundary from hit to L* (edge-by-edge); append vertices and L* to polyline
    static void walkToLstar(const amp::Obstacle2D& O, size_t startEdge,
                            const Eigen::Vector2d& hitPoint, Turn turn,
                            const Eigen::Vector2d& Lstar, size_t L_edge,
                            std::vector<Eigen::Vector2d>& polyline);

    // Append straight segment (dedups start point)
    static void appendSeg(std::vector<Eigen::Vector2d>& out,
                          const Eigen::Vector2d& a, const Eigen::Vector2d& b);

    Algo algo_;
    Turn turn_;
};
