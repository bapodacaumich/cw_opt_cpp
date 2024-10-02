#include "utils.hpp"

#include <casadi/casadi.hpp>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

casadi::MX cw_v_init(const std::vector<float>& start, casadi::MX& dx, casadi::MX& dy, casadi::MX& dz, casadi::MX& t){
    const float n=1.1288e-3;

    const float& x = start[0];
    const float& y = start[1];
    const float& z = start[2];

    casadi::MX cnt = cos(n * t);
    casadi::MX snt = sin(n * t);
    casadi::MX sigma3 = cnt * cnt;
    casadi::MX sigma2 = snt * snt;
    casadi::MX sigma1 = 4.0f * sigma3 - 8.0f * cnt + 4.0f * sigma2 - 3.0f * n * t * snt + 4.0f;

    casadi::MX vx = -n * (2.0f * dy - 2.0f * y - 2.0f * dy * cnt - 4.0f * dx * snt + 2.0f * y * cnt + 4.0f * x * snt + 3.0f * dx * n * t - 3.0f * n * t * x * cnt)/sigma1;
    casadi::MX vy = -n * (8.0f * x - 2.0f * dx + 2.0f * dx * cnt - dy * snt - 14.0f * x * cnt + y * snt + 6.0f * x * sigma3 + 6.0f * x * sigma2 - 6.0f * n * t * x * snt)/sigma1;
    casadi::MX vz = n * (dz - z * cnt) / snt;

    std::vector<casadi::MX> v_init_vec = {vx, vy, vz};
    casadi::MX v_init = casadi::MX::vertcat(v_init_vec);
    return v_init;
}

casadi::MX cw_v_init(casadi::MX& x, casadi::MX& y, casadi::MX& z, const std::vector<float>& end, casadi::MX& t){
    const float n=1.1288e-3;

    const float& dx = end[0];
    const float& dy = end[1];
    const float& dz = end[2];

    casadi::MX cnt = cos(n * t);
    casadi::MX snt = sin(n * t);
    casadi::MX sigma3 = cnt * cnt;
    casadi::MX sigma2 = snt * snt;
    casadi::MX sigma1 = 4.0f * sigma3 - 8.0f * cnt + 4.0f * sigma2 - 3.0f * n * t * snt + 4.0f;

    casadi::MX vx = -n * (2.0f * dy - 2.0f * y - 2.0f * dy * cnt - 4.0f * dx * snt + 2.0f * y * cnt + 4.0f * x * snt + 3.0f * dx * n * t - 3.0f * n * t * x * cnt)/sigma1;
    casadi::MX vy = -n * (8.0f * x - 2.0f * dx + 2.0f * dx * cnt - dy * snt - 14.0f * x * cnt + y * snt + 6.0f * x * sigma3 + 6.0f * x * sigma2 - 6.0f * n * t * x * snt)/sigma1;
    casadi::MX vz = n * (dz - z * cnt) / snt;

    std::vector<casadi::MX> v_init_vec = {vx, vy, vz};
    casadi::MX v_init = casadi::MX::vertcat(v_init_vec);
    return v_init;
}


casadi::MX cw_v_init(const std::vector<float>& start, const std::vector<float>& end, casadi::MX& t){
    const float n=1.1288e-3;

    const float& dx = end[0];
    const float& dy = end[1];
    const float& dz = end[2];
    const float& x = start[0];
    const float& y = start[1];
    const float& z = start[2];

    casadi::MX cnt = cos(n * t);
    casadi::MX snt = sin(n * t);
    casadi::MX sigma3 = cnt * cnt;
    casadi::MX sigma2 = snt * snt;
    casadi::MX sigma1 = 4.0f * sigma3 - 8.0f * cnt + 4.0f * sigma2 - 3.0f * n * t * snt + 4.0f;

    casadi::MX vx = -n * (2.0f * dy - 2.0f * y - 2.0f * dy * cnt - 4.0f * dx * snt + 2.0f * y * cnt + 4.0f * x * snt + 3.0f * dx * n * t - 3.0f * n * t * x * cnt)/sigma1;
    casadi::MX vy = -n * (8.0f * x - 2.0f * dx + 2.0f * dx * cnt - dy * snt - 14.0f * x * cnt + y * snt + 6.0f * x * sigma3 + 6.0f * x * sigma2 - 6.0f * n * t * x * snt)/sigma1;
    casadi::MX vz = n * (dz - z * cnt) / snt;

    std::vector<casadi::MX> v_init_vec = {vx, vy, vz};
    casadi::MX v_init = casadi::MX::vertcat(v_init_vec);
    return v_init;
}

casadi::DM cw_v_init(const std::vector<float>& start, const std::vector<float>& end, casadi::DM& t){
    const float n=1.1288e-3;

    const float& dx = end[0];
    const float& dy = end[1];
    const float& dz = end[2];
    const float& x = start[0];
    const float& y = start[1];
    const float& z = start[2];

    casadi::DM cnt = cos(n * t);
    casadi::DM snt = sin(n * t);
    casadi::DM sigma3 = cnt * cnt;
    casadi::DM sigma2 = snt * snt;
    casadi::DM sigma1 = 4.0f * sigma3 - 8.0f * cnt + 4.0f * sigma2 - 3.0f * n * t * snt + 4.0f;

    casadi::DM vx = -n * (2.0f * dy - 2.0f * y - 2.0f * dy * cnt - 4.0f * dx * snt + 2.0f * y * cnt + 4.0f * x * snt + 3.0f * dx * n * t - 3.0f * n * t * x * cnt)/sigma1;
    casadi::DM vy = -n * (8.0f * x - 2.0f * dx + 2.0f * dx * cnt - dy * snt - 14.0f * x * cnt + y * snt + 6.0f * x * sigma3 + 6.0f * x * sigma2 - 6.0f * n * t * x * snt)/sigma1;
    casadi::DM vz = n * (dz - z * cnt) / snt;

    std::vector<casadi::DM> v_init_vec = {vx, vy, vz};
    casadi::DM v_init = casadi::DM::vertcat(v_init_vec);
    return v_init;
}

casadi::MX cw_v_end(const std::vector<float>& start, const std::vector<float>& end, casadi::MX& t){
    const float n=1.1288e-3;

    const float& dx = end[0];
    const float& dy = end[1];
    const float& dz = end[2];
    const float& x = start[0];
    const float& z = start[2];

    casadi::MX cnt = cos(n * t);
    casadi::MX snt = sin(n * t);

    casadi::MX vx = dx * cnt + 2.0f * dy * snt + 3.0f * n * x * snt;
    casadi::MX vy = dy * (4.0f * cnt - 3.0f) - 2.0f * dx * snt + 6.0f * n * x * (cnt - 1.0f);
    casadi::MX vz = dz * cnt - n * z * snt;

    std::vector<casadi::MX> v_end_vec = {vx, vy, vz};
    casadi::MX v_end = casadi::MX::vertcat(v_end_vec);
    return v_end;
}

std::vector<casadi::MX> cw_pose(const std::vector<float>& start, casadi::MX& v_start, casadi::MX& t){
    // get pose along clohessy wiltshire trajectory
    const float& n=1.1288e-3;

    const float& x = start[0];
    const float& y = start[1];
    const float& z = start[2];

    casadi::MX vx = v_start(0);
    casadi::MX vy = v_start(1);
    casadi::MX vz = v_start(2);

    casadi::MX cnt = cos(n * t);
    casadi::MX snt = sin(n * t);

    casadi::MX xe = vx * snt / n - x * (3 * cnt - 4) - 2 * vy * (cnt - 1) / n;
    casadi::MX ye = y + x * (6 * snt - 6 * n * t) + vy * (4 * snt - 3 * n * t) / n + 2 * vx * (cnt - 1) / n;
    casadi::MX ze = z * cnt + vz * snt / n;

    std::vector<casadi::MX> pose_vec = {xe, ye, ze};
    return pose_vec;
}

std::vector<casadi::DM> cw_pose(const std::vector<float>& start, casadi::DM& v_start, casadi::DM& t){
    // get pose along clohessy wiltshire trajectory
    const float& n=1.1288e-3;

    const float& x = start[0];
    const float& y = start[1];
    const float& z = start[2];

    casadi::DM vx = v_start(0);
    casadi::DM vy = v_start(1);
    casadi::DM vz = v_start(2);

    casadi::DM cnt = cos(n * t);
    casadi::DM snt = sin(n * t);

    casadi::DM xe = vx * snt / n - x * (3 * cnt - 4) - 2 * vy * (cnt - 1) / n;
    casadi::DM ye = y + x * (6 * snt - 6 * n * t) + vy * (4 * snt - 3 * n * t) / n + 2 * vx * (cnt - 1) / n;
    casadi::DM ze = z * cnt + vz * snt / n;

    std::vector<casadi::DM> pose_vec = {xe, ye, ze};
    return pose_vec;
}

std::vector<casadi::MX> cw_pose(casadi::MX& x, casadi::MX& y, casadi::MX& z, casadi::MX& v_start, casadi::MX& t){
    // get pose along clohessy wiltshire trajectory
    const float& n=1.1288e-3;

    casadi::MX vx = v_start(0);
    casadi::MX vy = v_start(1);
    casadi::MX vz = v_start(2);

    casadi::MX cnt = cos(n * t);
    casadi::MX snt = sin(n * t);

    casadi::MX xe = vx * snt / n - x * (3 * cnt - 4) - 2 * vy * (cnt - 1) / n;
    casadi::MX ye = y + x * (6 * snt - 6 * n * t) + vy * (4 * snt - 3 * n * t) / n + 2 * vx * (cnt - 1) / n;
    casadi::MX ze = z * cnt + vz * snt / n;

    std::vector<casadi::MX> pose_vec = {xe, ye, ze};
    return pose_vec;
}

casadi::MX compute_path_cost(casadi::MX& T, const std::vector<std::vector<float>>& knot_points, bool square){
    /*
    compute the path cost for a given trajectory
    args:
    - T: casadi::MX, time vector
    - knot_points: std::vector<std::vector<float>>, list of knot points
    - square: bool, whether to square the cost
    */

    // Save velocity at end of every segment for dv computation
    casadi::MX v_last = casadi::MX::zeros(3, 1);

    // running dv sum (magnitude)
    casadi::MX dv_tot = casadi::MX::zeros(1, 1);
    for (size_t i = 0; i < knot_points.size(); i++)
    {
        casadi::MX dv = casadi::MX::zeros(1, 1);
        if (i == knot_points.size()-1) // last iteration, add stopping velocity
        {
            for (size_t j = 0; j < 3; j++)
            {
                dv += v_last(j) * v_last(j);
            }
        }
        else
        {
            casadi::MX t = T(i,0);
            casadi::MX v_init = cw_v_init(knot_points[i], knot_points[i+1], t);
            for (size_t j = 0; j < 3; j++)
            {
                dv += (v_init(j) - v_last(j)) * (v_init(j) - v_last(j));
            }
            if (!square)
            {
                dv = sqrt(dv);
            }
            v_last = cw_v_end(knot_points[i], knot_points[i+1], t);
        }
        dv_tot += dv;
    }
    return dv_tot;
}

bool loadCSV(const std::string& filename, std::vector<std::vector<float>>& data, int rowlen){
    /*
    load a csv file into a vector of vectors
    args:
    - filename: std::string, path to csv file
    - data: std::vector<std::vector<float>>, output data
    */
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line))
    {
        std::vector<float> row;
        std::stringstream ss(line);
        std::string cell;
        int num_cells = 0;
        while (std::getline(ss, cell, ',') && num_cells < rowlen)
        {
            try {
                row.push_back(std::stof(cell)); // Convert string to float and add to row
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid number format: " << cell << std::endl;
                return false;
            } catch (const std::out_of_range& e) {
                std::cerr << "Number out of range: " << cell << std::endl;
                return false;
            }
            ++num_cells;
        }
        data.push_back(row);
    }

    file.close();
    return true;
}

void saveCSV(const std::string& filename, const std::vector<std::vector<float>>& data) {
    /*
    * Save 2d std::vector float to a csv file
    * @param filename: std::string, path to save file
    * @param data: std::vector<std::vector<float>>, data to save
    */


    // Open the file in output mode
    std::ofstream file(filename);
    
    // Check if the file was opened successfully
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    // Iterate over the rows
    size_t n_rows = data.size();
    for (size_t i = 0; i < n_rows; ++i) {
        // Iterate over the columns
        size_t n_cols = data[i].size();
        for (size_t j = 0; j < n_cols; ++j) {
            file << std::fixed << std::setprecision(6) << data[i][j]; // Write value
            if (j < n_cols - 1) {
                file << ","; // Separate values with a comma
            }
        }
        file << "\n"; // End of row
    }

    // Close the file
    file.close();
}


void saveCSV(const std::string& filename, const casadi::DM& data) {
    /*
    * Save casadi::DM to a csv file
    * @param filename: std::string, path to save file
    * @param data: casadi::DM, data to save
    */


    // Open the file in output mode
    std::ofstream file(filename);
    
    // Check if the file was opened successfully
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    
    // Iterate over the rows
    size_t n_rows = data.size1();
    for (size_t i = 0; i < n_rows; ++i) {
        // Iterate over the columns
        size_t n_cols = data.size2();
        for (size_t j = 0; j < n_cols; ++j) {
            file << std::fixed << std::setprecision(16) << data(i,j); // Write value
            if (j < n_cols - 1) {
                file << ","; // Separate values with a comma
            }
        }
        file << "\n"; // End of row
    }

    // Close the file
    file.close();
}

std::string firstFileWithPrefix(const std::string& directory, const std::string& prefix) {
    /*
    * Find the first file in a directory with a given prefix
    * @param directory: std::string, path to directory
    * @param prefix: std::string, prefix to search for
    */
    namespace fs = std::filesystem;
    for (const auto& entry : std::filesystem::directory_iterator(directory)) {
        if (entry.is_regular_file() && entry.path().filename().string().find(prefix) == 0) {
            return entry.path().string();
        }
    }
    return "";
}

void print2DVector(const std::vector<std::vector<float>>& vec) {
    for (const auto& row : vec) {
        for (float value : row) {
            std::cout << std::fixed << std::setprecision(2) << value << ' '; // Print with fixed precision
        }
        std::cout << '\n'; // New line after each row
    }
}

void displayProgressBar(double progress, int width) {
    // Clear the current line
    std::cout << "\r";

    // Calculate the number of '#' characters
    int pos = static_cast<int>(width * progress);
    
    // Draw the progress bar
    std::cout << "[";
    for (int i = 0; i < width; ++i) {
        if (i < pos) 
            std::cout << "#";
        else 
            std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(1) << (progress * 100) << "%";
    
    // Flush the output to ensure it updates immediately
    std::cout.flush();
}

casadi::DM get_intermediate_points_init(const std::vector<std::vector<float>>& knot_points, const casadi::DM& T_init){
    // get intermediate points for warm start
    casadi::DM X_init = casadi::DM::zeros(knot_points.size()-1, 3);
    for (size_t i = 0; i < knot_points.size()-1; i++)
    {
        casadi::DM t = T_init(2 * i,0) + T_init(2 * i + 1,0);
        casadi::DM v_init = cw_v_init(knot_points[i], knot_points[i+1], t);
        casadi::DM t0 = T_init(2 * i,0);
        std::vector<casadi::DM> ip = cw_pose(knot_points[i], v_init, t0);
        X_init(i,0) = ip[0];
        X_init(i,1) = ip[1];
        X_init(i,2) = ip[2];
    }
    return X_init;
}