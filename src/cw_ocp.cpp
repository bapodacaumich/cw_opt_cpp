#include "utils.hpp"

#include <casadi/casadi.hpp>
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
    std::cout << "Locality: " << argv[2] << ", int repr: " << locality << std::endl;
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

    // load file contents into knot_points
    std::vector<std::vector<float>> knot_points;
    loadCSV(matching_file, knot_points);

    print2DVector(knot_points);

    size_t n_knots = knot_points.size();

    auto opti = casadi::Opti();

    // Constraints
    casadi::MX T = opti.variable(n_knots-1, 1);
    opti.subject_to(sum1(sum2(T)) < tmax);

    for (size_t i = 0; i < n_knots-1; i++)
    {
        opti.subject_to(T(i,0) > 0);
    }

    // Objective
    casadi::MX dv_tot = compute_path_cost(T, knot_points, true);
    opti.minimize(dv_tot);

    // Warm start
    casadi::DM T_init = casadi::DM::ones(n_knots-1, 1);
    for (size_t i = 0; i < n_knots-1; i++)
    {
        T_init(i,0) = tmax / (n_knots-1) / 2;
    }
    opti.set_initial(T, T_init);

    casadi::Dict p_opts;
    casadi::Dict s_opts;
    s_opts["print_level"] = 7;
    s_opts["tol"] = 1e-3;
    s_opts["max_iter"] = 10000;

    opti.solver("ipopt", p_opts, s_opts);

    casadi::OptiSol sol = opti.solve();

    std::cout << sol.value(T) << std::endl; // print solution

    return 0;
}