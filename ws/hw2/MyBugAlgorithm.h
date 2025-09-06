#pragma once

#include "AMPCore.h"
#include "hw/HW2.h"
#include <vector>
#include <Eigen/Dense>

class MyBugAlgorithm : public amp::BugAlgorithm {
public:
    enum Turn { LEFT, RIGHT };

    // Default is a right-hand-on-wall traversal
    explicit MyBugAlgorithm(Turn turn = RIGHT) : turn_(turn) {}

    // The main planning method
    virtual amp::Path2D plan(const amp::Problem2D& problem) override;

private:
    // Struct to hold information about a collision
    struct Hit {
        bool hit = false;
        size_t obsIdx = 0;
        size_t edgeIdx = 0;
        Eigen::Vector2d point;
    };

    // Helper Functions
    static std::vector<Eigen::Vector2d> vertsOf(const amp::Obstacle2D& obs);
    static bool segSegIntersect(const Eigen::Vector2d& p, const Eigen::Vector2d& q,
                                const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                                double& t, double& u, Eigen::Vector2d& x);
    static inline double dgoal(const Eigen::Vector2d& p, const Eigen::Vector2d& g) {
        return (p - g).norm();
    }
    static Eigen::Vector2d closestPointOnSegment(const Eigen::Vector2d& A,
                                                 const Eigen::Vector2d& B,
                                                 const Eigen::Vector2d& P);
    static Hit findFirstHit(const Eigen::Vector2d& s, const Eigen::Vector2d& g,
                            const std::vector<amp::Obstacle2D>& obstacles);
    
    // Surveys the entire boundary to find L*
    static void surveyBoundaryForLstar(const amp::Obstacle2D& O, Turn turn,
                                       const Eigen::Vector2d& qg,
                                       Eigen::Vector2d& Lstar, size_t& L_edgeIdx);
    
    // Generates the path from the hit point to the leave point (L*)
    static void travelToLeavePoint(const amp::Obstacle2D& O, size_t hit_edgeIdx,
                                   const Eigen::Vector2d& hit_point, Turn turn,
                                   const Eigen::Vector2d& Lstar, size_t L_edgeIdx,
                                   std::vector<Eigen::Vector2d>& polyline);

    static Eigen::Vector2d getOutwardNormal(const std::vector<Eigen::Vector2d>& V,
                                            size_t edgeIdx, const Eigen::Vector2d& point);
    static void appendSegment(std::vector<Eigen::Vector2d>& out,
                              const Eigen::Vector2d& a, const Eigen::Vector2d& b);

    Turn turn_;
};
