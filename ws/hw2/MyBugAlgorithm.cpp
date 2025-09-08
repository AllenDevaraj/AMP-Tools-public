#include "MyBugAlgorithm.h"
#include <algorithm>
#include <limits>
#define _USE_MATH_DEFINES
#include <math.h>

// =================================================================
// ALGORITHM PARAMETERS - You can tune these!
// =================================================================
static constexpr double STEP_SIZE = 0.1;
static constexpr double WALL_FOLLOW_DISTANCE = 0.2;
static constexpr double SENSOR_RANGE = 0.5;
static constexpr double GOAL_RADIUS = 0.2;
static constexpr double HIT_POINT_RADIUS = 0.15; // How close to get to the start to count as a full loop


amp::Path2D MyBugAlgorithm::plan(const amp::Problem2D& problem) {
    amp::Path2D path;
    path.waypoints.push_back(problem.q_init);
    
    Eigen::Vector2d current_pos = problem.q_init;
    int loop_guard = 0;

    while ((current_pos - problem.q_goal).norm() > GOAL_RADIUS && loop_guard++ < 2000) {
        // STATE 1: MOVE TO GOAL
        Eigen::Vector2d heading_to_goal = (problem.q_goal - current_pos).normalized();
        SensorHit front_hit = checkSensor(current_pos, heading_to_goal, (problem.q_goal - current_pos).norm(), problem.obstacles);

        if (!front_hit.found) {
            // Path to goal is clear, we are done
            path.waypoints.push_back(problem.q_goal);
            break;
        }
        
        // Add path to the hit point
        Eigen::Vector2d hit_point = current_pos + heading_to_goal * front_hit.distance;
        path.waypoints.push_back(hit_point);
        current_pos = hit_point;

        // STATE 2: WALL FOLLOW (BUG 1)
        // Turn left 90 degrees to start the wall follow
        Eigen::Rotation2Dd turn_left(M_PI / 2.0);
        Eigen::Vector2d initial_follow_heading = turn_left * heading_to_goal;

        // This single function call now handles the entire Bug 1 process:
        // discovering the boundary, mapping it, and finding L*.
        WallFollowResult result = performBug1WallFollow(current_pos, initial_follow_heading, problem);

        // Now, find L* in the path we just discovered
        size_t l_star_index = 0;
        double min_dist = std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < result.discovery_path.size(); ++i) {
            double dist = (result.discovery_path[i] - problem.q_goal).norm();
            if (dist < min_dist) {
                min_dist = dist;
                l_star_index = i;
            }
        }
        
        // Append the path from the hit point to L*
        for (size_t i = 0; i <= l_star_index; ++i) {
            path.waypoints.push_back(result.discovery_path[i]);
        }
        
        current_pos = result.discovery_path[l_star_index];

        // Leave the obstacle
        Eigen::Vector2d last_step = (result.discovery_path[l_star_index] - result.discovery_path[l_star_index - 1]).normalized();
        Eigen::Rotation2Dd turn_right(-M_PI / 2.0);
        Eigen::Vector2d outward_normal = turn_right * last_step; // Normal is perpendicular to the direction of travel
        current_pos += outward_normal * 0.1; // Nudge away from the wall
        path.waypoints.push_back(current_pos);
    }

    return path;
}


MyBugAlgorithm::WallFollowResult MyBugAlgorithm::performBug1WallFollow(const Eigen::Vector2d& hit_point, const Eigen::Vector2d& initial_heading, const amp::Problem2D& problem) {
    WallFollowResult result;
    result.closest_point_to_goal = hit_point;
    double min_dist_to_goal = (hit_point - problem.q_goal).norm();

    Eigen::Vector2d current_pos = hit_point;
    Eigen::Vector2d heading = initial_heading;

    int follow_steps = 0;
    // Loop until we make a full circle and get back to the hit_point
    while (follow_steps++ < 10000) {
        // Reactive steering logic
        Eigen::Rotation2Dd turn_right(-M_PI / 2.0);
        Eigen::Rotation2Dd turn_left(M_PI / 2.0);
        Eigen::Vector2d right_dir = turn_right * heading;
        Eigen::Vector2d front_dir = heading;

        SensorHit right_sensor = checkSensor(current_pos, right_dir, SENSOR_RANGE, problem.obstacles);
        SensorHit front_sensor = checkSensor(current_pos, front_dir, WALL_FOLLOW_DISTANCE, problem.obstacles);

        if (front_sensor.found) {
            heading = turn_left * heading; // Inner corner, turn left
        } else if (!right_sensor.found || right_sensor.distance > WALL_FOLLOW_DISTANCE * 1.5) {
            heading = turn_right * heading; // Outer corner, turn right
        } else {
            double error = WALL_FOLLOW_DISTANCE - right_sensor.distance;
            Eigen::Rotation2Dd correction(error * 0.5); // Corrective steering
            heading = correction * heading;
        }

        current_pos += STEP_SIZE * heading.normalized();
        result.discovery_path.push_back(current_pos);

        // While discovering, keep track of the closest point to the goal (L*)
        double current_dist_to_goal = (current_pos - problem.q_goal).norm();
        if (current_dist_to_goal < min_dist_to_goal) {
            min_dist_to_goal = current_dist_to_goal;
            result.closest_point_to_goal = current_pos;
        }

        // Check if we have completed a full lap
        if (follow_steps > 50 && (current_pos - hit_point).norm() < HIT_POINT_RADIUS) {
            break; // We are back at the start, survey is complete.
        }
    }
    
    return result;
}


// =================================================================
// SENSOR AND INTERSECTION HELPERS
// =================================================================

MyBugAlgorithm::SensorHit MyBugAlgorithm::checkSensor(const Eigen::Vector2d& pos, const Eigen::Vector2d& dir, double range,
                                                      const std::vector<amp::Obstacle2D>& obstacles) {
    SensorHit result;
    Eigen::Vector2d ray_end = pos + dir.normalized() * range;

    for (const auto& obs : obstacles) {
        auto vertices = obs.verticesCW();
        for (size_t i = 0; i < vertices.size(); ++i) {
            Eigen::Vector2d p1 = vertices[i];
            Eigen::Vector2d p2 = vertices[(i + 1) % vertices.size()];
            
            double t, u;
            if (segSegIntersect(pos, ray_end, p1, p2, t, u)) {
                if (t >= -1e-9 && t <= 1.0 + 1e-9) {
                    double dist = t * range;
                    if (dist < result.distance) {
                        result.distance = dist;
                        result.found = true;
                    }
                }
            }
        }
    }
    return result;
}

bool MyBugAlgorithm::segSegIntersect(const Eigen::Vector2d& p, const Eigen::Vector2d& q,
                                     const Eigen::Vector2d& a, const Eigen::Vector2d& b,
                                     double& t, double& u) {
    const Eigen::Vector2d r = q - p;
    const Eigen::Vector2d s = b - a;
    const double denom = r.x() * s.y() - r.y() * s.x();
    const double num_t = (a - p).x() * s.y() - (a - p).y() * s.x();
    const double num_u = (a - p).x() * r.y() - (a - p).y() * r.x();

    if (std::abs(denom) < 1e-9) return false;

    t = num_t / denom;
    u = num_u / denom;

    if (t >= -1e-9 && t <= 1.0 + 1e-9 && u >= -1e-9 && u <= 1.0 + 1e-9) {
        return true;
    }
    return false;
}