#pragma once

#include "AMPCore.h"
#include <Eigen/Dense>
#include <vector>

class MyBugAlgorithm {
public:
    // Bug type and Direction
    enum class BugMode    { Bug1, Bug2 };
    enum class FollowSide { Right, Left };

    explicit MyBugAlgorithm(BugMode mode = BugMode::Bug1,
                            FollowSide side = FollowSide::Right)
        : mode_(mode), side_(side) {}

    amp::Path2D plan(const amp::Problem2D& problem);

    // Constants
    static constexpr double STEP_SIZE       = 0.04;
    static constexpr double CONTACT_TOL     = 1e-3;
    static constexpr double GOAL_RADIUS     = 0.20;
    static constexpr double HIT_RADIUS      = 0.22;
    static constexpr double WALL_DIST       = 0.22;
    static constexpr double SENSOR_RANGE    = 2.00;
    static constexpr double IMPROVE_MIN     = 1e-4;
    static constexpr double LEAVE_EPS       = 1e-2;
    static constexpr double GOAL_EPS        = 3e-3;
    static constexpr double MIN_LOOP_ARC    = 4.0;

    static constexpr int    MAX_OUTER_ITERS = 1800;
    static constexpr int    MAX_STEPS  = 300000;

    struct SensorHit {
        bool   found    = false;
        double distance = 0.0;
    };

    struct WallFollowResult {
        Eigen::Vector2d              Lstar;  
        std::vector<Eigen::Vector2d> path_gone;      // Path travelled during Wall Follow
        double                       arc_length = 0.0;       // distance traced along boundary
    };

    // Optional: quick access to last traced boundary arc
    double lastBoundaryArc() const { return last_boundary_arc_; }

private:
    //  Modes
    BugMode    mode_ = BugMode::Bug1;
    FollowSide side_ = FollowSide::Right;

    SensorHit checkSensor(const Eigen::Vector2d& pos,
                          const Eigen::Vector2d& dir,
                          double range,
                          const std::vector<amp::Obstacle2D>& obstacles);

    // Bug-1 wall follow
    WallFollowResult Bug1WallFollow(const Eigen::Vector2d& hit_point,
                                    const Eigen::Vector2d& initial_heading_right,
                                    const amp::Problem2D& problem);

    // Plans
    amp::Path2D planBug1(const amp::Problem2D& problem);
    amp::Path2D planBug2(const amp::Problem2D& problem);

    struct HitInfo {
        Eigen::Vector2d H;       // hitpoint
        double          rho_H;   // ||H - goal||
        Eigen::Vector2d m_hat;   // unit direction of M-line (q_init -> q_goal)
        Eigen::Vector2d A;       // q_init
        Eigen::Vector2d B;       // q_goal
    };

    bool onMLineCloser(const Eigen::Vector2d& p,
                       const amp::Problem2D& P,
                       const HitInfo& hit) const;

    void Bug2WallFollow(Eigen::Vector2d& q,
                                    const HitInfo& hit,
                                    const amp::Problem2D& P,
                                    std::vector<Eigen::Vector2d>& out_trace);

    // bookkeeping for distance traced during the last boundary follow
    double last_boundary_arc_ = 0.0;
};
