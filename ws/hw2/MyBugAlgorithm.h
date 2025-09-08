#pragma once

#include "AMPCore.h"
#include "hw/HW2.h"
#include <vector>
#include <Eigen/Dense>

class MyBugAlgorithm : public amp::BugAlgorithm {
public:
    MyBugAlgorithm() {}

    virtual amp::Path2D plan(const amp::Problem2D& problem) override;

private:
    // Struct to hold sensor information
    struct SensorHit {
        bool found = false;
        double distance = std::numeric_limits<double>::infinity();
    };

    // Struct to return the results of the wall-following discovery process
    struct WallFollowResult {
        Eigen::Vector2d closest_point_to_goal; // This will be L*
        std::vector<Eigen::Vector2d> discovery_path; // The path taken during circumnavigation
    };

    // Sensor simulation function
    SensorHit checkSensor(const Eigen::Vector2d& pos, const Eigen::Vector2d& dir, double range,
                          const std::vector<amp::Obstacle2D>& obstacles);
    
    // Main function for the new reactive boundary following
    WallFollowResult performBug1WallFollow(const Eigen::Vector2d& hit_point, const Eigen::Vector2d& initial_heading,
                                           const amp::Problem2D& problem);

    // Basic intersection helper
    static bool segSegIntersect(const Eigen::Vector2d& p, const Eigen::Vector2d& q,
                                const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                                double& t, double& u);
};