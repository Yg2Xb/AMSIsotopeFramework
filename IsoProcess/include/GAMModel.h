#ifndef GAMMODEL_H
#define GAMMODEL_H

#include <vector>
#include "TNamed.h"
#include "TH2D.h"
#include <iostream>

class GAMModel : public TNamed {
    public:
        std::vector<double> buffer;      // Buffer to store the weights
        std::vector<std::vector<float>> centroids;
        std::vector<double> scaling;    
        std::vector<double> means;
        TH2D index_correction;

        GAMModel() : TNamed() {} // 使用 TNamed 的默认构造函数

        void retrieve_weights(float *args);
        void retrieve_weights(double *args);
        double eval(float *args);
        double eval(double *args);
        double get_index_correction(float x, float y);

        ClassDef(GAMModel, 1);
};

#endif