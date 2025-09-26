#include "AMPCore.h"
#include <Eigen/Dense>
#include <algorithm> 

double cross_product(const Eigen::Vector2d& O, const Eigen::Vector2d& A, const Eigen::Vector2d& B) {
    return (A.x() - O.x()) * (B.y() - O.y()) - (A.y() - O.y()) * (B.x() - O.x());
}

std::vector<Eigen::Vector2d> convexHull(std::vector<Eigen::Vector2d>& points) {
    int n = points.size();
    if (n <= 3) return points;
    std::sort(points.begin(), points.end(), [](const Eigen::Vector2d& a, const Eigen::Vector2d& b) {
        return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y());
    });
    std::vector<Eigen::Vector2d> hull;
    for (int i = 0; i < n; ++i) {
        while (hull.size() >= 2 && cross_product(hull[hull.size()-2], hull.back(), points[i]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(points[i]);
    }
    for (int i = n - 2, t = hull.size() + 1; i >= 0; i--) {
        while (hull.size() >= t && cross_product(hull[hull.size()-2], hull.back(), points[i]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(points[i]);
    }
    hull.pop_back();
    return hull;
}

// Helper function to perform the Minkowski sum
amp::Polygon minkowskiSum(const amp::Polygon& polygonA, const amp::Polygon& polygonB) {
    std::vector<Eigen::Vector2d> sum_vertices;
    for (const auto& vA : polygonA.verticesCCW()) {
        for (const auto& vB : polygonB.verticesCCW()) {
            sum_vertices.push_back(vA + vB);
        }
    }
    std::vector<Eigen::Vector2d> hull_vertices = convexHull(sum_vertices);
    return amp::Polygon(hull_vertices);
}

// Main function to drive the visualization
int main(int argc, char** argv) {
    std::vector<Eigen::Vector2d> obstacle_vertices = {Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.0, 2.0), Eigen::Vector2d(0.0, 2.0)};
    amp::Polygon workspace_obstacle(obstacle_vertices);

    std::vector<Eigen::Vector2d> robot_vertices = {Eigen::Vector2d(0.0, 0.0), Eigen::Vector2d(1.0, 2.0), Eigen::Vector2d(0.0, 2.0)};
    amp::Polygon robot_shape(robot_vertices);

    int num_slices = 12;
    std::vector<amp::Polygon> cspace_slices;
    std::vector<double> angles;

    for (int i = 0; i < num_slices; ++i) {
        double theta = i * (2.0 * M_PI / num_slices);
        angles.push_back(theta);

        std::vector<Eigen::Vector2d> rotated_vertices;
        Eigen::Matrix2d rotation_matrix;
        rotation_matrix << cos(theta), -sin(theta), sin(theta), cos(theta);
        for (const auto& vertex : robot_shape.verticesCCW()) rotated_vertices.push_back(rotation_matrix * vertex);
        amp::Polygon rotated_robot(rotated_vertices);

        std::vector<Eigen::Vector2d> reflected_vertices;
        for (const auto& vertex : rotated_robot.verticesCCW()) reflected_vertices.push_back(-vertex);
        amp::Polygon reflected_robot(reflected_vertices);

        amp::Polygon cspace_obstacle_slice = minkowskiSum(workspace_obstacle, reflected_robot);
        cspace_slices.push_back(cspace_obstacle_slice);
    }

    amp::Visualizer::makeFigure(cspace_slices, angles);
    // amp::Visualizer::saveFigures(); 

    return 0;
}
