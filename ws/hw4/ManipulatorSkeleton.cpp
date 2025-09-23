#include "ManipulatorSkeleton.h"
#include <cmath> 

MyManipulator2D::MyManipulator2D()
    : LinkManipulator2D({1.0, 1.0})
{}

// Override this method for implementing forward kinematics
Eigen::Vector2d MyManipulator2D::getJointLocation(const amp::ManipulatorState& state, uint32_t joint_index) const {
    // Defensive check to handle improperly sized states
    const amp::ManipulatorState* state_ptr = &state;
    amp::ManipulatorState temp_state;
    if (state.size() != nLinks()) {
        temp_state.resize(nLinks());
        temp_state.setZero();
        state_ptr = &temp_state;
    }

    // BUG FIX: The previous loop logic was incorrect. This is a simpler and correct implementation.
    Eigen::Vector2d joint_location(0.0, 0.0);
    double cumulative_angle = 0.0;
    
    // To find the location of joint `j`, we sum the vectors of links 0 through j-1.
    for (uint32_t i = 0; i < joint_index; ++i) {
        cumulative_angle += (*state_ptr)(i);
        joint_location.x() += m_link_lengths[i] * cos(cumulative_angle);
        joint_location.y() += m_link_lengths[i] * sin(cumulative_angle);
    }

    return joint_location;
}

// Override this method for implementing inverse kinematics
amp::ManipulatorState MyManipulator2D::getConfigurationFromIK(const Eigen::Vector2d& end_effector_location) const {
    // This implementation is for a 2-link manipulator
    if (nLinks() != 2) {
        amp::ManipulatorState joint_angles(nLinks());
        joint_angles.setZero();
        return joint_angles;
    }

    amp::ManipulatorState state(2);
    const double x = end_effector_location.x();
    const double y = end_effector_location.y();
    const double l1 = m_link_lengths[0];
    const double l2 = m_link_lengths[1];

    // Check if the target is reachable
    const double distance_squared = x*x + y*y;
    if (distance_squared > (l1 + l2) * (l1 + l2) || distance_squared < (l1 - l2) * (l1 - l2)) {
        state.setZero();
        return state;
    }

    // Calculate theta2 using the Law of Cosines (elbow up)
    const double cos_theta2 = (x*x + y*y - l1*l1 - l2*l2) / (2.0 * l1 * l2);
    const double theta2 = -acos(cos_theta2);

    // Calculate theta1
    const double theta1 = atan2(y, x) - atan2(l2 * sin(theta2), l1 + l2 * cos(theta2));

    state << theta1, theta2;
    return state;
}