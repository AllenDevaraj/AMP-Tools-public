#pragma once

#include "AMPCore.h"
#include <Eigen/Dense>
#include <vector>

class MyBugAlgorithm {
public:
    // ===== Bug modes =====
    enum class BugMode { Bug1, Bug2 };

    explicit MyBugAlgorithm(BugMode mode = BugMode::Bug1) : mode_(mode) {}

    // Planner API expected by HW2::grade / Visualizer
    amp::Path2D plan(const amp::Problem2D& problem);

    // ---------- Tunables ----------
    static constexpr double STEP_SIZE       = 0.04;
    static constexpr double CONTACT_TOL     = 1e-3;
    static constexpr double GOAL_RADIUS     = 0.20;
    static constexpr double HIT_RADIUS      = 0.22;
    static constexpr double WALL_DIST       = 0.22;
    static constexpr double SENSOR_RANGE    = 2.00;
    static constexpr double IMPROVE_MIN     = 1e-4;
    static constexpr double LEAVE_EPS       = 1e-2;
    static constexpr double GOAL_BIAS_EPS   = 3e-3;
    static constexpr double MIN_LOOP_ARC    = 4.0;

    static constexpr int    MAX_OUTER_ITERS = 1800;
    static constexpr int    MAX_WALK_STEPS  = 300000;

    // ---------- PODs ----------
    struct SensorHit {
        bool   found    = false;
        double distance = 0.0;
    };

    struct WallFollowResult {
        Eigen::Vector2d              closest_point_to_goal;  // L*
        std::vector<Eigen::Vector2d> discovery_path;         // boundary samples H … ~H
    };

private:
    // ===== Mode =====
    BugMode mode_ = BugMode::Bug1;

    // ===== Shared helpers (implemented in .cpp) =====
    SensorHit checkSensor(const Eigen::Vector2d& pos,
                          const Eigen::Vector2d& dir,
                          double range,
                          const std::vector<amp::Obstacle2D>& obstacles);

    // Bug-1 wall follow
    WallFollowResult Bug1WallFollow(const Eigen::Vector2d& hit_point,
                                    const Eigen::Vector2d& initial_heading_right,
                                    const amp::Problem2D& problem);

    // Split plans
    amp::Path2D planBug1(const amp::Problem2D& problem);
    amp::Path2D planBug2(const amp::Problem2D& problem);

    // ===== Bug-2 helpers =====
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

    void followRightUntilMLineLeave(Eigen::Vector2d& q,
                                    const HitInfo& hit,
                                    const amp::Problem2D& P,
                                    std::vector<Eigen::Vector2d>& out_trace);
};
