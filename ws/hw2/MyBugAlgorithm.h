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
    static constexpr double STEP_SIZE          = 0.04;   // was 0.06
    static constexpr double CONTACT_TOL        = 1e-3;
    static constexpr double GOAL_RADIUS        = 0.20;
    static constexpr double HIT_RADIUS         = 0.22;   // was 0.15 (easier close near hit)
    static constexpr double WALL_DIST          = 0.22;   // was 0.30 (better for narrow gaps)
    static constexpr double SENSOR_RANGE       = 2.00;   // was 1.00 (more reliable re-acquire)
    static constexpr double IMPROVE_MIN        = 1e-4;   // was 1e-3 (don’t declare blocked too early)
    static constexpr double LEAVE_EPS          = 1e-2;   // was 3e-3 (peel off a tad more)
    static constexpr double GOAL_BIAS_EPS      = 3e-3;
    static constexpr double MIN_LOOP_ARC       = 4.0;    // was 10.0 (random obstacles can be small)
    static constexpr int    MAX_OUTER_ITERS    = 1800;
    static constexpr int    MAX_WALK_STEPS     = 300000; // tiny headroom


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

// ====== BUG-2: separate class (does not modify your Bug-1) =================