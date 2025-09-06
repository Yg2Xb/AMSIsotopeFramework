#include "GAMModel.h"
#include <cmath>

ClassImp(GAMModel);

void GAMModel::retrieve_weights(float *args) {
    int dimensions = scaling.size();
    buffer.clear();
    buffer.push_back(1); // This is the constant part of the fit

    for(int i = 0; i < dimensions; i++) 
        buffer.push_back(args[i]); // Linear part

    for(int i = 0; i < dimensions; i++) 
        for(auto &v: centroids[i]) 
            buffer.push_back(pow(fabs(v-args[i])/scaling[i], 3));
}

void GAMModel::retrieve_weights(double *args) {
    int dimensions = scaling.size();
    buffer.clear();
    buffer.push_back(1); // This is the constant part of the fit

    for(int i = 0; i < dimensions; i++) 
        buffer.push_back(args[i]); // Linear part

    for(int i = 0; i < dimensions; i++) 
        for(auto &v: centroids[i]) 
            buffer.push_back(pow(fabs(v-args[i])/scaling[i], 3));
}

double GAMModel::eval(float *args) {
    retrieve_weights(args);
    double v = 0;
    for(unsigned int i = 0; i < means.size(); i++) 
        v += buffer[i]*means[i];
    return v;
}

double GAMModel::eval(double *args) {
    retrieve_weights(args);
    double v = 0;
    for(unsigned int i = 0; i < means.size(); i++) 
        v += buffer[i]*means[i];
    return v;
}

double GAMModel::get_index_correction(float x, float y) {
    if(index_correction.GetNbinsX()==1 || index_correction.GetNbinsY()==1) 
        return 1;
    int bx = index_correction.GetXaxis()->FindBin(x);
    int by = index_correction.GetYaxis()->FindBin(y);
    return 1-index_correction.GetBinContent(bx,by);
}