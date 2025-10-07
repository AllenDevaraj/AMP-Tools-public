#pragma once

#include "AMPCore.h"
#include "hw/HW6.h"      // base classes (GridCSpace2D, *CSConstructor, WaveFrontAlgorithm)
#include "HelpfulClass.h"

#include <utility>
#include <cstddef>
#include <Eigen/Dense>

// ------------------------------------------------------------
// Custom grid C-space that stores its own dims/bounds and
// exposes helpers used by the wave-front implementation.
// ------------------------------------------------------------
class MyGridCSpace2D : public amp::GridCSpace2D {
public:
    MyGridCSpace2D(std::size_t nx, std::size_t ny,
                   double x_min, double x_max,
                   double y_min, double y_max)
    : amp::GridCSpace2D(nx, ny, x_min, x_max, y_min, y_max),
      m_nx(nx), m_ny(ny),
      m_xmin(x_min), m_xmax(x_max),
      m_ymin(y_min), m_ymax(y_max) {}

    // Map a workspace point to a grid cell (ix,iy)
    std::pair<std::size_t, std::size_t>
    getCellFromPoint(double x0, double x1) const override;

    // Centers of cells in world coordinates
    Eigen::Vector2d cellCenter(std::size_t ix, std::size_t iy) const {
        const double dx = (m_xmax - m_xmin) / static_cast<double>(m_nx);
        const double dy = (m_ymax - m_ymin) / static_cast<double>(m_ny);
        const double cx = m_xmin + (static_cast<double>(ix) + 0.5) * dx;
        const double cy = m_ymin + (static_cast<double>(iy) + 0.5) * dy;
        return Eigen::Vector2d(cx, cy);
    }

    // Dimensions (number of cells)
    std::size_t width()  const { return m_nx; }
    std::size_t height() const { return m_ny; }

    // Bounds (if needed elsewhere)
    double xmin() const { return m_xmin; }
    double xmax() const { return m_xmax; }
    double ymin() const { return m_ymin; }
    double ymax() const { return m_ymax; }

private:
    std::size_t m_nx, m_ny;
    double m_xmin, m_xmax, m_ymin, m_ymax;
};

// ------------------------------------------------------------
// Point-agent C-space constructor (Exercise 1 only)
// ------------------------------------------------------------
class MyPointAgentCSConstructor : public amp::PointAgentCSConstructor {
public:
    explicit MyPointAgentCSConstructor(std::size_t cells_per_dim)
    : m_cells_per_dim(cells_per_dim) {}

    std::unique_ptr<amp::GridCSpace2D>
    construct(const amp::Environment2D& env) override;

private:
    std::size_t m_cells_per_dim;
};

// ------------------------------------------------------------
// Manipulator constructor (NOT used in Ex1; minimal stub)
// ------------------------------------------------------------
class MyManipulatorCSConstructor : public amp::ManipulatorCSConstructor {
public:
    explicit MyManipulatorCSConstructor(std::size_t cells_per_dim)
    : m_cells_per_dim(cells_per_dim) {}

    std::unique_ptr<amp::GridCSpace2D>
    construct(const amp::LinkManipulator2D& manipulator,
              const amp::Environment2D& env) override;

private:
    std::size_t m_cells_per_dim;
};

// ------------------------------------------------------------
// Concrete WaveFront (since base is abstract in your build)
//  - 4-connected BFS on grid centers
//  - endpoints forced to q_init, q_goal
// ------------------------------------------------------------
class MyWaveFrontAlgorithm : public amp::WaveFrontAlgorithm {
public:
    amp::Path2D planInCSpace(const Eigen::Vector2d& q_init,
                             const Eigen::Vector2d& q_goal,
                             const amp::GridCSpace2D& grid_cspace,
                             bool isManipulator) override;
};
