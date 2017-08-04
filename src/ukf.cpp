#include "ukf.h"
#include "Eigen/Dense"
#include "tools.h"
#include <iostream>

#define EPS 0.001 // Just a small number

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    // initial state vector
    x_ = VectorXd(5);

    // initial covariance matrix
    P_ = MatrixXd(5, 5);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = .5;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = .5;

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    previous_timestamp_ = 0;

    // State dimension
    n_x_ = x_.size();

    // Augmented state dimension
    n_aug_ = n_x_ + 2;

    // Number of sigma points
    n_sig_ = 2 * n_aug_ + 1;

    // Set the predicted sigma points matrix dimentions
    Xsig_pred_ = MatrixXd(n_x_, n_sig_);

    // Sigma point spreading parameter
    lambda_ = 3 - n_aug_;

    n_z_ = 3;

    // Weights of sigma points
    weights_ = VectorXd(n_sig_);

    // Measurement noise covariance matrices initialization
    R_radar_ = MatrixXd(3, 3);
    R_radar_ << std_radr_ * std_radr_, 0, 0,
            0, std_radphi_ * std_radphi_, 0,
            0, 0, std_radrd_ * std_radrd_;

    R_lidar_ = MatrixXd(2, 2);
    R_lidar_ << std_laspx_ * std_laspx_, 0,
            0, std_laspy_ * std_laspy_;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {

        //state covariance matrix P
        P_ <<
           1, 0, 0, 0, 0,
                0, 1, 0, 0, 0,
                0, 0, 1, 0, 0,
                0, 0, 0, 1, 0,
                0, 0, 0, 0, 1;

        // first measurement
        cout << "UKF: " << endl;
        x_ << 1, 1, 1, 1, 1;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */

            float rho = meas_package.raw_measurements_[0]; // range
            float phi = meas_package.raw_measurements_[1]; // bearing
            float rho_dot = meas_package.raw_measurements_[2]; // velocity of rho

            // Coordinates conversion from polar to cartesian
            float px = rho * cos(phi);
            float py = rho * sin(phi);
            float vx = rho_dot * cos(phi);
            float vy = rho_dot * sin(phi);
            float v = sqrt(vx * vx + vy * vy);

            x_ << px, py, v, 0, 0;

        } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
            // We don't know velocities from the first measurement of the LIDAR, so, we use zeros
            x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

            // Deal with the special case initialisation problems
            if (fabs(x_(0)) < EPS and fabs(x_(1)) < EPS) {
                x_(0) = EPS;
                x_(1) = EPS;
            }
        }

        // Initialize weights
        weights_(0) = lambda_ / (lambda_ + n_aug_);
        for (int i = 1; i < weights_.size(); i++) {
            weights_(i) = 0.5 / (n_aug_ + lambda_);
        }

        previous_timestamp_ = meas_package.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    //compute the time elapsed between the current and previous measurements
    double dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
    previous_timestamp_ = meas_package.timestamp_;

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    Prediction(dt);

    /*****************************************************************************
     *  Update
     ****************************************************************************/
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
        UpdateRadar(meas_package);
    }

    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
        UpdateLidar(meas_package);
    }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
    //calculate square root of P
    MatrixXd A = P_.llt().matrixL();

    //set first column of sigma point matrix
    Xsig_pred_.col(0) = x_;

    //set remaining sigma points
    for (int i = 0; i < n_x_; i++) {
        Xsig_pred_.col(i + 1) = x_ + sqrt(lambda_ + n_x_) * A.col(i);
        Xsig_pred_.col(i + 1 + n_x_) = x_ - sqrt(lambda_ + n_x_) * A.col(i);
    }

    //create augmented mean vector
    VectorXd x_aug = VectorXd(7);

    //create augmented state covariance
    MatrixXd P_aug = MatrixXd(7, 7);

    //create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

    x_aug.head(5) = x_;
    x_aug[5] = 0;
    x_aug[6] = 0;

    P_aug.fill(0.0);
    P_aug.topLeftCorner(5, 5) = P_;

    P_aug(5, 5) = std_a_ * std_a_;
    P_aug(6, 6) = std_yawdd_ * std_yawdd_;

    MatrixXd L = P_aug.llt().matrixL();

    //set first column of sigma point matrix
    Xsig_aug.col(0) = x_aug;

    //set remaining sigma points
    for (int i = 0; i < n_aug_; i++) {
        Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
        Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
    }

    for (int i = 0; i < 2 * n_aug_ + 1; i++) {
        //extract values for better readability
        double p_x = Xsig_aug(0, i);
        double p_y = Xsig_aug(1, i);
        double v = Xsig_aug(2, i);
        double yaw = Xsig_aug(3, i);
        double yawd = Xsig_aug(4, i);
        double nu_a = Xsig_aug(5, i);
        double nu_yawdd = Xsig_aug(6, i);

        //predicted state values
        double px_p, py_p;

        //avoid division by zero
        if (fabs(yawd) > EPS) {
            px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
            py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
        } else {
            px_p = p_x + v * delta_t * cos(yaw);
            py_p = p_y + v * delta_t * sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd * delta_t;
        double yawd_p = yawd;

        //add noise
        px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
        py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
        v_p = v_p + nu_a * delta_t;

        yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
        yawd_p = yawd_p + nu_yawdd * delta_t;

        //write predicted sigma point into right column
        Xsig_pred_(0, i) = px_p;
        Xsig_pred_(1, i) = py_p;
        Xsig_pred_(2, i) = v_p;
        Xsig_pred_(3, i) = yaw_p;
        Xsig_pred_(4, i) = yawd_p;

        //predicted state mean
        x_.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
            x_ = x_ + weights_(i) * Xsig_pred_.col(i);
        }

        //predict state covariance matrix
        P_.fill(0.0);
        for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

            // state difference
            VectorXd x_diff = Xsig_pred_.col(i) - x_;

            //angle normalization
            while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
            while (x_diff(3) < -M_PI) x_diff(3) += 2. * M_PI;

            P_ = P_ + weights_(i) * x_diff * x_diff.transpose();

        }
    }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    // Set measurement dimension
    int n_z = 2;

    // Create matrix for sigma points in measurement space
    // Transform sigma points into measurement space
    MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_sig_);

    UpdateUKF(meas_package, Zsig, n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    MatrixXd Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);

    //transform sigma points into measurement space
    for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

        // extract values for better readibility
        double p_x = Xsig_pred_(0, i);
        double p_y = Xsig_pred_(1, i);
        double v = Xsig_pred_(2, i);
        double yaw = Xsig_pred_(3, i);

        double v1 = cos(yaw) * v;
        double v2 = sin(yaw) * v;

        // measurement model
        Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);                           //r
        Zsig(1, i) = atan2(p_y, p_x);                                       //phi
        Zsig(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y);   //r_dot
    }

    UpdateUKF(meas_package, Zsig, n_z_);
}

/**
 *  Angle normalization to [-Pi, Pi]
 */
void UKF::Norma(double *ang) {
    while (*ang > M_PI) *ang -= 2. * M_PI;
    while (*ang < -M_PI) *ang += 2. * M_PI;
}

void UKF::UpdateUKF(MeasurementPackage meas_package, MatrixXd Zsig, int n_z) {
    // Mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);
    z_pred = Zsig * weights_;

    //measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z, n_z);
    S.fill(0.0);

    for (int i = 0; i < n_sig_; i++) {
        // Residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // Angle normalization
        Norma(&(z_diff(1)));

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }
    // Add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z, n_z);

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) { // Radar
        R = R_radar_;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) { // Lidar
        R = R_lidar_;
    }

    S = S + R;

    // Create matrix for cross correlation Tc
    MatrixXd Tc = MatrixXd(n_x_, n_z);

    // Calculate cross correlation matrix
    Tc.fill(0.0);

    for (int i = 0; i < n_sig_; i++) {
        //residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        if (meas_package.sensor_type_ == MeasurementPackage::RADAR) { // Radar
            // Angle normalization
            Norma(&(z_diff(1)));
        }

        // State difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;

        // Angle normalization
        Norma(&(x_diff(3)));

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // Measurements
    VectorXd z = meas_package.raw_measurements_;

    //Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    // Residual
    VectorXd z_diff = z - z_pred;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) { // Radar
        // Angle normalization
        Norma(&(z_diff(1)));
    }

    // Update state mean and covariance matrix
    x_ = x_ + K * z_diff;
    P_ = P_ - K * S * K.transpose();

    VectorXd temp = meas_package.raw_measurements_ - z_pred;
    nis = temp.transpose() * S.inverse() * temp;

    // Update NIS graph
    tools.updateGraph(nis, meas_package.sensor_type_);
}