#pragma once

// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW4.h"

// Derive the amp::GridCSpace2D class and override the missing method
class MyGridCSpace2D : public amp::GridCSpace2D {
    public:
        // MODIFIED: The constructor now initializes its own member variables in addition to calling the base class.
        MyGridCSpace2D(std::size_t x0_cells, std::size_t x1_cells, double x0_min, double x0_max, double x1_min, double x1_max)
            : amp::GridCSpace2D(x0_cells, x1_cells, x0_min, x0_max, x1_min, x1_max), // Call base class constructor
              m_x0_min(x0_min), m_x0_max(x0_max), m_x1_min(x1_min), m_x1_max(x1_max),
              m_x0_cells(x0_cells), m_x1_cells(x1_cells)
        {}

        // Override this method for determining which cell a continuous point belongs to
        virtual std::pair<std::size_t, std::size_t> getCellFromPoint(double x0, double x1) const override;

        // NEW: Public getter methods so other classes (like MyManipulatorCSConstructor) can access the grid properties.
        double getX0_min() const { return m_x0_min; }
        double getX0_max() const { return m_x0_max; }
        double getX1_min() const { return m_x1_min; }
        double getX1_max() const { return m_x1_max; }
        std::size_t getX0_cells() const { return m_x0_cells; }
        std::size_t getX1_cells() const { return m_x1_cells; }

    private:
        // NEW: Private member variables to store the grid properties.
        double m_x0_min, m_x0_max, m_x1_min, m_x1_max;
        std::size_t m_x0_cells, m_x1_cells;
};

// Derive the HW4 ManipulatorCSConstructor class and override the missing method
class MyManipulatorCSConstructor : public amp::ManipulatorCSConstructor {
    public:
        // To make things easy, add the number of cells as a ctor param so you can easily play around with it
        MyManipulatorCSConstructor(std::size_t cells_per_dim) : m_cells_per_dim(cells_per_dim) {}

        // Override this method for computing all of the boolean collision values for each cell in the cspace
        virtual std::unique_ptr<amp::GridCSpace2D> construct(const amp::LinkManipulator2D& manipulator, const amp::Environment2D& env) override;

    private:
        std::size_t m_cells_per_dim;
};