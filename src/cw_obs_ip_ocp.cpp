#include "utils.hpp"
#include "constraints.hpp"

#include <casadi/casadi.hpp>
#include <chrono>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

int main(int argc, char* argv[]){
    // argc: number of arguments
    // argv: array of arguments

    std::cout << "Number of args: " << argc << std::endl;

    char* vgd = NULL;
    bool locality = true;
    float tmax = 10.0f;

    if (argc > 1) {
        vgd = argv[1];

        if (argc > 2) {
            locality = (argv[2][0] == 't');

            if (argc > 3) {
                tmax = std::stof(argv[3]);
                std::stringstream ss(argv[3]);
                ss >> tmax;

                if (ss.fail()) {
                    std::cerr << "T max parsing failed: " << argv[3] << std::endl;
                    return 0;
                }
            }
        }
    }

    std::cout << "VGD: " << vgd << std::endl;
    std::cout << "Locality: " << argv[2] << std::endl;
    std::cout << "T max: " << tmax << std::endl;

    std::string locality_str = "";
    if (locality) {locality_str = "_local";}

    // find file that matches vgd, locality, and knot directory
    std::string knotdir = "../data/knot_points/";
    std::string matching_file = firstFileWithPrefix(knotdir, vgd+locality_str);
    if (matching_file == "") {
        std::cerr << "No matching file found in directory: " << knotdir << std::endl;
        return 0;
    }

    // load in station normals and centroids
    std::vector<std::vector<std::vector<float>>> station_normals;
    std::vector<std::vector<std::vector<float>>> station_centroids;

    size_t num_obstacles = 15;
    for (size_t i = 0; i < num_obstacles; i++){
        std::string normal_file = "../data/model/normal_centroid/obs" + std::to_string(i+1) + "_normals.txt";
        std::string centroid_file = "../data/model/normal_centroid/obs" + std::to_string(i+1) + "_points.txt";

        std::vector<std::vector<float>> normals;
        std::vector<std::vector<float>> centroids;

        loadCSV(normal_file, normals);
        loadCSV(centroid_file, centroids);

        station_normals.push_back(normals);
        station_centroids.push_back(centroids);
    }

    // load file contents into knot_points
    std::vector<std::vector<float>> knot_points;
    loadCSV(matching_file, knot_points);

    size_t n_knots = knot_points.size();

    auto opti = casadi::Opti();

    // Variables
    size_t n_drift = (n_knots -1) * 2;
    casadi::MX T = opti.variable(n_drift, 1);
    casadi::MX X = opti.variable(n_knots - 1, 3);

    // Constraints
    opti.subject_to(sum1(sum2(T)) < tmax);

    for (size_t i = 0; i < n_drift; i++)
    {
        opti.subject_to(T(i,0) > 0);
    }

    auto start = std::chrono::high_resolution_clock::now();

    enforce_station_keepout(station_normals, station_centroids, opti, knot_points, X, T);

    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;

    std::cout << "Time to enforce station keepout constraints: " << elapsed.count() << " seconds " << std::endl;

    // Objective
    casadi::MX dv_tot = compute_path_cost(T, knot_points, true);
    opti.minimize(dv_tot);

    // Warm start
    casadi::DM T_init = casadi::DM::ones(n_drift, 1);
    for (size_t i = 0; i < n_drift; i++)
    {
        T_init(i,0) = tmax / n_drift / 2;
    }
    opti.set_initial(T, T_init);

    casadi::DM X_init = get_intermediate_points_init(knot_points, T_init);
    opti.set_initial(X, X_init);

    casadi::Dict p_opts;
    casadi::Dict s_opts;
    s_opts["print_level"] = 7;
    s_opts["tol"] = 1e-3;
    s_opts["max_iter"] = 10000;

    std::cout << "Setting up OCP..." << std::endl;
    opti.solver("ipopt", p_opts, s_opts);

    casadi::OptiSol sol = opti.solve();

    std::cout << sol.value(T) << std::endl; // print solution
    std::cout << sol.value(X) << std::endl; // print solution

    saveCSV("../data/solution/obs_ip/" + std::string(vgd) + locality_str + "_t.csv", sol.value(T));
    saveCSV("../data/solution/obs_ip/" + std::string(vgd) + locality_str + "_x.csv", sol.value(X));

    return 0;
}