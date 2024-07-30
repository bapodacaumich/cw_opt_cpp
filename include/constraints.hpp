#ifndef CONSTRAINTS_HPP
#define CONSTRAINTS_HPP

#include <casadi/casadi.hpp>
#include <vector>

void enforce_convex_hull_from_points(
    const std::vector<std::vector<float>>& normals,
    const std::vector<std::vector<float>>& points,
    casadi::Opti& opti,
    std::vector<std::vector<casadi::MX>>& X,
    float min_station_distance
    );

void enforce_station_keepout(
    const std::vector<std::vector<std::vector<float>>>& station_normals,
    const std::vector<std::vector<std::vector<float>>>& station_points, 
    casadi::Opti& opti, 
    const std::vector<std::vector<float>>& knot_points, 
    casadi::MX& T, 
    float min_station_distance=0.0f, 
    size_t nT=3
    );

void enforce_station_keepout(
    const std::vector<std::vector<std::vector<float>>>& station_normals,
    const std::vector<std::vector<std::vector<float>>>& station_points, 
    casadi::Opti& opti, 
    const std::vector<std::vector<float>>& knot_points, 
    casadi::MX& X,
    casadi::MX& T, 
    float min_station_distance=0.0f, 
    size_t nT=3
    );

#endif // CONSTRAINTS_HPP