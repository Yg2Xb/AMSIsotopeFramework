#ifndef MODELMANAGER_H
#define MODELMANAGER_H

#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TSystem.h>
#include "GAMModel.h" 

enum Rad
{
    AGL = 0,
    NAF = 1
};

class ModelManager
{
public:
    static GAMModel model[2][2];
    static ModelManager *head;

    // Initialization
    static ModelManager *init(TString filename_data, TString filename_mc);
    static float corrected_beta(float beta,
                              Rad radiator,
                              float run,
                              float charge,
                              float x_rad,
                              float y_rad,
                              float theta_rad,
                              float phi_rad,
                              float nreflected,
                              float nhits,
                              int mc = 0);
};

#endif // MODELMANAGER_H