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

    // ---------- Tunables ----------
    static constexpr double STEP_SIZE          = 0.06;   // motion step (m-line & wall-follow)
    static constexpr double CONTACT_TOL        = 1e-3;   // bisection tolerance for contact
    static constexpr double GOAL_RADIUS        = 0.20;   // done if within this of goal
    static constexpr double HIT_RADIUS         = 0.15;   // loop closure proximity to hit
    static constexpr double WALL_DIST          = 0.30;   // desired right-hand stand-off from wall
    static constexpr double SENSOR_RANGE       = 1.00;   // ray range for side sensor
    static constexpr double IMPROVE_MIN        = 1e-3;   // L* must be this much closer than hit
    static constexpr double LEAVE_EPS          = 3e-3;   // outward peel-off nudge
    static constexpr double GOAL_BIAS_EPS      = 3e-3;   // tiny push toward goal on leave
    static constexpr double MIN_LOOP_ARC       = 10.0;   // must traverse at least this much boundary
    static constexpr int    MAX_OUTER_ITERS    = 1800;
    static constexpr int    MAX_WALK_STEPS     = 250000;

    // ---------- Sensor output ----------
    struct SensorHit {
        bool   found    = false;                                      // hit along ray?
        double distance = std::numeric_limits<double>::infinity();    // from ray start to first contact
    };

    // ---------- Circumnavigation result ----------
    struct WallFollowResult {
        Vec2                      closest_point_to_goal; // L*
        std::vector<Vec2>         discovery_path;        // recorded loop (hit → ... → near hit)
    };

    // ---------- Boolean “range” sensor ----------
    // March in STEP_SIZE; when next step would collide, bisect to CONTACT_TOL.
    SensorHit checkSensor(const Vec2& pos, const Vec2& dir, double range,
                          const std::vector<amp::Obstacle2D>& obstacles);

    // ---------- Full Bug-1 survey with RIGHT-hand follow (clockwise) ----------
    WallFollowResult performBug1WallFollow(const Vec2& hit_point,
                                           const Vec2& initial_heading_right,
                                           const amp::Problem2D& problem);
};
