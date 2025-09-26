#pragma once

#include "AMPCore.h"
#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <algorithm>

#define _USE_MATH_DEFINES
#include <math.h>

namespace MotionPlanningHelpers {
    using Vec = Eigen::Vector2d;

    inline double dist(const Vec& a, const Vec& b) {
        return (a - b).norm();
    }

    inline Vec unit(const Vec& v) {
        double n = v.norm();
        if (n > 0.0) return v / n;
        return Vec(1.0, 0.0);
    }

    inline double angle(const Vec& v) {
        return std::atan2(v.y(), v.x());
    }

    inline Vec ang2vec(double th) {
        return Vec(std::cos(th), std::sin(th));
    }

    inline double angDiff(double a, double b) {
        double d = std::fmod(a - b + M_PI, 2.0 * M_PI);
        if (d < 0) d += 2.0 * M_PI;
        return d - M_PI;
    }

    inline void append(std::vector<Vec>& path, const Vec& p, double tol = 1e-9) {
        if (path.empty() || (path.back() - p).norm() > tol) {
            path.push_back(p);
        }
    }

    class CollisionChecker {
    private:
        static constexpr double EPS = 1e-9;

        static std::vector<Vec> getClockwiseVertices(const amp::Obstacle2D& obstacle) {
            std::vector<Vec> V = obstacle.verticesCW();
            if (V.size() < 3) return V;

            double A = 0.0;
            for (size_t i = 0; i < V.size(); ++i) {
                const auto& p = V[i];
                const auto& q_cur = V[(i + 1) % V.size()];
                A += p.x() * q_cur.y() - q_cur.x() * p.y();
            }
            if (A > 0) std::reverse(V.begin(), V.end());
            return V;
        }

        static bool isPointInsidePolygon(const Vec& point, const std::vector<Vec>& polygon_vertices) {
            bool inside = false;
            const size_t n = polygon_vertices.size();
            for (size_t i = 0, j = n - 1; i < n; j = i++) {
                const Vec& A = polygon_vertices[j];
                const Vec& B = polygon_vertices[i];
                const bool cond = ((A.y() > point.y()) != (B.y() > point.y())) &&
                                  (point.x() < (B.x() - A.x()) * (point.y() - A.y()) / (B.y() - A.y()) + A.x());
                if (cond) inside = !inside;
            }
            return inside;
        }
        
        static bool checkSegmentIntersection(const Vec& p1, const Vec& q1, const Vec& p2, const Vec& q2) {
            const Vec r = q1 - p1;
            const Vec s = q2 - p2;
            const double denom = r.x() * s.y() - r.y() * s.x();
            const double num_t = (p2.x() - p1.x()) * s.y() - (p2.y() - p1.y()) * s.x();
            const double num_u = (p2.x() - p1.x()) * r.y() - (p2.y() - p1.y()) * r.x();

            if (std::abs(denom) < EPS) return false;

            const double t = num_t / denom;
            const double u = num_u / denom;

            return (t >= -EPS && t <= 1.0 + EPS && u >= -EPS && u <= 1.0 + EPS);
        }

    public:
        /**
         * @brief Checks if a point is in collision with any obstacle.
         * @param point The point to check.
         * @param obstacles A vector of obstacles in the environment.
         * @return True if the point is inside any obstacle, false otherwise.
         */
        static bool isPointInCollision(const Vec& point, const std::vector<amp::Obstacle2D>& obstacles) {
            for (const auto& obs : obstacles) {
                auto vertices = getClockwiseVertices(obs);
                if (vertices.size() < 3) continue;
                if (isPointInsidePolygon(point, vertices)) return true; // In collision
            }
            return false; // Not in collision (free)
        }

        /**
         * @brief Checks if a line segment is in collision with any obstacle.
         * @param p1 The start point of the segment.
         * @param p2 The end point of the segment.
         * @param obstacles A vector of obstacles in the environment.
         * @return True if the segment crosses any obstacle edge or has an endpoint inside an obstacle, false otherwise.
         */
        static bool isSegmentInCollision(const Vec& p1, const Vec& p2, const std::vector<amp::Obstacle2D>& obstacles) {
            if (isPointInCollision(p1, obstacles) || isPointInCollision(p2, obstacles)) {
                return true;
            }

            for (const auto& obs : obstacles) {
                auto vertices = getClockwiseVertices(obs);
                const size_t n = vertices.size();
                if (n < 2) continue;
                for (size_t i = 0; i < n; ++i) {
                    if (checkSegmentIntersection(p1, p2, vertices[i], vertices[(i + 1) % n])) {
                        return true;
                    }
                }
            }
            return false;
        }
    };

} 
