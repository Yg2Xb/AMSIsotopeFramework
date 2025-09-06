#include "ModelManager.h"
#include <cmath>

// Initialize static members
GAMModel ModelManager::model[2][2];
ModelManager* ModelManager::head = nullptr;

ModelManager* ModelManager::init(TString filename_data, TString filename_mc) {
    if (head != nullptr) {
        std::cerr << "ModelManager::init -- already called." << std::endl;
        return head;
    }

    ///////////////// DATA CORRECTION /////////////////
    {
        TFile fin(filename_data);
        GAMModel *m = (GAMModel *)fin.Get("model_z_agl");
        if (!m) {
            std::cerr << "ModelManager::init -- File " << filename_data << " does not contains model_z_agl" << std::endl;
            exit(1);
        }
        model[0][0] = *m;

        m = (GAMModel *)fin.Get("model_z_naf");
        if (!m) {
            std::cerr << "ModelManager::init -- File " << filename_data << " does not contains model_z_naf" << std::endl;
            exit(1);
        }
        model[0][1] = *m;
    }

    ///////////////// MC CORRECTION /////////////////
    {
        TFile fin(filename_mc);
        GAMModel *m = (GAMModel *)fin.Get("model_z_agl");
        if (!m) {
            std::cerr << "ModelManager::init -- File " << filename_mc << " does not contains model_z_agl" << std::endl;
            exit(1);
        }
        model[1][0] = *m;

        m = (GAMModel *)fin.Get("model_z_naf");
        if (!m) {
            std::cerr << "ModelManager::init -- File " << filename_mc << " does not contains model_z_naf" << std::endl;
            exit(1);
        }
        model[1][1] = *m;
    }

    head = new ModelManager;
    return head;
}

float ModelManager::corrected_beta(float beta,
                                 Rad radiator,
                                 float run,
                                 float charge,
                                 float x_rad,
                                 float y_rad,
                                 float theta_rad,
                                 float phi_rad,
                                 float nreflected,
                                 float nhits,
                                 int mc) {
    float vx = sin(theta_rad) * cos(phi_rad);
    float vy = sin(theta_rad) * sin(phi_rad);

    float args[] = {
        nreflected / nhits,
        static_cast<float>(run - 1.30845e+09),
        charge
    };

    // SELECTION: THE CORRECTION DOES NOT APPLY FOR THESE
    if (std::min(fabs(x_rad), fabs(y_rad)) >= 40.5)
        return 0;
    if (x_rad * x_rad + y_rad * y_rad > 58.5 * 58.5)
        return 0;
    if (std::max(fabs(x_rad), fabs(y_rad)) > 28.5 && std::max(fabs(x_rad), fabs(y_rad)) < 29.5)
        return 0;

    if (mc == 1)
        args[1] = 0; // Turn off run number of MC
    double k = model[mc][radiator].eval(args);

    if (k == 0)
        std::cout << "WHAT!!!!!" << " " << args[0] << " " << args[1] << " " << args[2] << std::endl;

    if (!k)
        return 0;
    return beta * model[mc][radiator].get_index_correction(x_rad, y_rad) / k;
}