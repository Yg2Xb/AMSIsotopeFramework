#include <cmath>
#include "ModelManager.h"
#include "GAMModel.h"

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
        GAMModel *m_agl = (GAMModel *)fin.Get("model_z_agl");
        if (!m_agl) {
            std::cerr << "ModelManager::init -- File " << filename_data << " does not contains model_z_agl" << std::endl;
            exit(1);
        }
        // 使用 Clone() 方法进行深拷贝
        model[0][0] = *(GAMModel*)m_agl->Clone();

        GAMModel *m_naf = (GAMModel *)fin.Get("model_z_naf");
        if (!m_naf) {
            std::cerr << "ModelManager::init -- File " << filename_data << " does not contains model_z_naf" << std::endl;
            exit(1);
        }
        // 使用 Clone() 方法进行深拷贝
        model[0][1] = *(GAMModel*)m_naf->Clone();
    }

    ///////////////// MC CORRECTION /////////////////
    {
        TFile fin(filename_mc);
        GAMModel *m_agl = (GAMModel *)fin.Get("model_z_agl");
        if (!m_agl) {
            std::cerr << "ModelManager::init -- File " << filename_mc << " does not contains model_z_agl" << std::endl;
            exit(1);
        }
        // 使用 Clone() 方法进行深拷贝
        model[1][0] = *(GAMModel*)m_agl->Clone();

        GAMModel *m_naf = (GAMModel *)fin.Get("model_z_naf");
        if (!m_naf) {
            std::cerr << "ModelManager::init -- File " << filename_mc << " does not contains model_z_naf" << std::endl;
            exit(1);
        }
        // 使用 Clone() 进行深拷贝
        model[1][1] = *(GAMModel*)m_naf->Clone();
    }

    head = new ModelManager;
    return head;
}

void ModelManager::cleanup() {
    // 静态对象在程序退出时会自动清理
    // 这里只需清理动态分配的head
    if (head) {
        delete head;
        head = nullptr;
    }
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
                                    

    // Avoid overcorrecting outside of charge and time validity values  2025.12.15
    if(charge>9) charge=9;
    if(run>1.58013e+09) run=1.58013e+09;

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
