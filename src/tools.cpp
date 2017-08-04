#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {

    //reset the files and destroy the contents
    myfile = fopen("NIS_lidar", "w");
    fclose(myfile);

    myfile = fopen("NIS_radar", "w");
    fclose(myfile);

    gp_ = popen("gnuplot -persist", "w");

}

void Tools::updateGraph(double value, MeasurementPackage::SensorType type) {
    switch (type) {
        case MeasurementPackage::SensorType::LASER:
            myfile = fopen("NIS_lidar", "a");
            break;
        case MeasurementPackage::SensorType::RADAR:
            myfile = fopen("NIS_radar", "a");
            break;
    }

    if (myfile != NULL) {
        string toWrite = to_string(value) + '\n';
        fwrite(toWrite.c_str(), 1, toWrite.length(), myfile);
        fclose(myfile);
    }

    if (gp_ != NULL) {
        fprintf(gp_, "call '../plotScript.gp'\n");
        fflush(gp_);
    }
}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // Check the validity of the following inputs:
    // The estimation vector size should not be zero
    if (estimations.size() == 0) {
        cout << "Input is empty" << endl;
        return rmse;
    }

    // The estimation vector size should equal ground truth vector size
    if (estimations.size() != ground_truth.size()) {
        cout << "Invalid estimation or ground_truth. Data should have the same size" << endl;
        return rmse;
    }

    // Accumulate squared residuals
    for (unsigned int i = 0; i < estimations.size(); ++i) {
        VectorXd residual = estimations[i] - ground_truth[i];
        // Coefficient-wise multiplication
        residual = residual.array() * residual.array();
        rmse += residual;
    }

    // Calculate the mean
    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}
