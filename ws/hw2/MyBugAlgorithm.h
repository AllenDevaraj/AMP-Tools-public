#pragma once

#include "AMPCore.h"
#include "hw/HW2.h"
#include <Eigen/Dense>
#include <vector>
#include <limits>

class MyBugAlgorithm : public amp::BugAlgorithm {
public:
    MyBugAlgorithm() = default;

    // Main entry
    amp::Path2D plan(const amp::Problem2D& problem) override;

private:
    using Vec2 = Eigen::Vector2d;

    // ---------- Tunables (adjust if needed) ----------
    static constexpr double STEP_SIZE          = 0.10;   // marching step for m-line & wall follow
    static constexpr double CONTACT_TOL        = 1e-3;   // bisection tolerance to localize contact
    static constexpr double GOAL_RADIUS        = 0.20;   // done if inside this radius of goal
    static constexpr double HIT_RADIUS         = 0.15;   // loop closure threshold near hit point
    static constexpr double WALL_DIST          = 0.30;   // desired stand-off from wall (right side)
    static constexpr double SENSOR_RANGE       = 1.00;   // ray range for wall sensors
    static constexpr double IMPROVE_MIN        = 1e-3;   // L* must be at least this much closer than hit
    static constexpr double LEAVE_EPS          = 3e-3;   // outward peel-off nudge
    static constexpr double GOAL_BIAS_EPS      = 3e-3;   // tiny push toward goal on leave
    static constexpr int    MAX_OUTER_ITERS    = 1500;
    static constexpr int    MAX_WALK_STEPS     = 200000;

    // ---------- Sensor output ----------
    struct SensorHit {
        bool   found   = false;                                       // hit along ray?
        double distance = std::numeric_limits<double>::infinity();    // from ray start to first contact
    };

    // ---------- Circumnavigation result ----------
    struct WallFollowResult {
        Eigen::Vector2d              closest_point_to_goal; // L*
        std::vector<Eigen::Vector2d> discovery_path;        // polyline of the full loop (hit → ... → near hit)
    };

    // ---------- Boolean “range” sensor (planner never uses vertices directly) ----------
    // March in STEP_SIZE increments; when the *next* step would collide, bisect [a,b]
    // to the boundary with CONTACT_TOL accuracy and return that distance.
    SensorHit checkSensor(const Vec2& pos, const Vec2& dir, double range,
                          const std::vector<amp::Obstacle2D>& obstacles);

    // ---------- Full Bug-1 survey with RIGHT-hand wall follow (clockwise) ----------
    // Starting at hit_point with initial heading that keeps the wall on the RIGHT,
    // perform one full loop, recording the path and the best (closest-to-goal) L*.
    WallFollowResult performBug1WallFollow(const Vec2& hit_point,
                                           const Vec2& initial_heading_right,
                                           const amp::Problem2D& problem);
};
