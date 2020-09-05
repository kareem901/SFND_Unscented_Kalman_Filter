#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  std_a_ = 2.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  is_initialized_=false;
// State dimension 
  n_x_=5;
//Augmented State dimension
  n_aug_=7;
    // Define spreading Parameter
  lambda_ =3 -n_aug_;
  // predicted sigma points matrix
  Xsig_pred_ =MatrixXd(n_x_,2*n_aug_+1);

  //Vector for weights 
  weights_=VectorXd(2*n_aug_+1);
   //Start time 
  time_us_=0;
  // Current NIS for Radar 
  NIS_radar =0.0;

    // Current NIS for Laser 
  NIS_Laser =0.0;
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_)
  {
    if(meas_package.sensor_type_==MeasurementPackage::LASER)
    {
      float x=meas_package.raw_measurements_(0);
      float y=meas_package.raw_measurements_(1);
      x_<<x,y,0,0,0;
      P_<<std_laspx_*std_laspx_,0,0,0,0,
      0,std_laspy_*std_laspy_,0,0,0,
      0,0,1,0,0,
      0,0,0,1,0,
      0,0,0,0,1 ;
    }
    else
    {
      float rho=meas_package.raw_measurements_(0);
      float phi=meas_package.raw_measurements_(1);
      float rho_dot=meas_package.raw_measurements_(2);

      float x=rho*sin(phi);
      float y=rho*cos(phi);
      x_<<x,y,rho_dot,phi,0;
      P_<<std_radr_*std_radr_,0,0,0,0,
        0,std_radr_*std_radr_,0,0,0,
        0,0,std_radrd_*std_radrd_,0,0,
        0,0,0,std_radphi_*std_radphi_,0,
        0,0,0,0,1;

    }
    is_initialized_=true;
    time_us_=meas_package.timestamp_;
    return;
  }
  float dt=(meas_package.timestamp_-time_us_)/1000000.0;
  time_us_=meas_package.timestamp_;
  Prediction(dt);
   if(meas_package.sensor_type_==MeasurementPackage::LASER)
    {
       UpdateLidar(meas_package);
    }
    else
    {
       UpdateRadar(meas_package);
    }
    

}

void UKF::Prediction(double delta_t) {
//Generate Augmated Segma points 

//crate augmented mean vector
VectorXd x_aug=VectorXd(7);

//crate augmented state covariance  
MatrixXd P_aug=MatrixXd(7,7);

//crate Segma point matrix
MatrixXd Xsig_aug=MatrixXd(n_aug_,2*n_aug_+1);

//crate Augmented Mean state 
x_aug.head(5)=x_;
x_aug(5)=0;
x_aug(6)=0;

//crate Augmented Covariance Matrix
P_aug.fill(0.0);
P_aug.topLeftCorner(5,5)=P_;
P_aug(5,5)=std_a_*std_a_;
P_aug(6,6)=std_yawdd_*std_yawdd_;

//Create Square Root matrix
MatrixXd L=P_aug.llt().matrixL();

//Create Augmented Sigma Points
Xsig_aug.fill(0.0);
Xsig_aug.col(0)=x_aug;
for (int i=0;i<n_aug_;i++)
{
  Xsig_aug.col(i+1)=x_aug+sqrt(lambda_+n_aug_)*L.col(i);
  Xsig_aug.col(i+1+n_aug_)=x_aug-sqrt(lambda_+n_aug_)*L.col(i);
}


//Predict Segma Points
for (int i=0;i<2*n_aug_+1;i++)
{
  double p_x=Xsig_aug(0,i);
  double p_y=Xsig_aug(1,i);
  double v=Xsig_aug(2,i);
  double yaw=Xsig_aug(3,i);
  double yawd=Xsig_aug(4,i);
  double nu_a=Xsig_aug(5,i);
  double nu_yawd=Xsig_aug(6,i);

  //Predict State Values 
  double px_p ,py_p;

  //Avoid Division by Zero
  if (fabs(yawd)>0.001)
  {
    px_p=p_x+v/yawd *(sin(yaw+yawd*delta_t)-sin(yaw));
    py_p=p_y +v/yawd*(cos(yaw)-cos(yaw+yawd*delta_t));
  }
  else
  {
    px_p=p_x+v*delta_t*cos(yaw);
     py_p=p_y+v*delta_t*sin(yaw);
  }

  double v_p=v;
  double yaw_p =yaw +yawd*delta_t;
  double yawd_p =yawd;

    //add noise 
    px_p=px_p +0.5*nu_a*delta_t*delta_t*cos(yaw);
     py_p=py_p +0.5*nu_a*delta_t*delta_t*sin(yaw);
     v_p=v_p+nu_a*delta_t;

     yaw_p=yaw_p+0.5*nu_yawd*delta_t*delta_t;
     yawd_p=yawd_p+nu_yawd*delta_t;

     //Write Prediction Segma Point
     Xsig_pred_(0,i)=px_p;
     Xsig_pred_(1,i)=py_p;
     Xsig_pred_(2,i)=v_p;
     Xsig_pred_(3,i)=yaw_p;
     Xsig_pred_(4,i)=yawd_p;

}
//Predicated Mean and Covariance 

//Set Weights 
double weight_0=lambda_/(lambda_+n_aug_);
weights_(0)=weight_0;
for (int i=1;i<2*n_aug_+1;i++)
{
  double weight=0.5/(n_aug_+lambda_);
  weights_(i)=weight;
}
//Predicted State mean 
x_.fill(0.0);
for (int i=0;i<2*n_aug_+1;i++)
{
  x_=x_+weights_(i)*Xsig_pred_.col(i);
}

//Predicted State Covariance matrix
P_.fill(0.0);
for (int i=0;i<2*n_aug_+1;i++)
{
  VectorXd x_diff=Xsig_pred_.col(i)-x_;
    while (x_diff(3) > M_PI) x_diff(3)-=2*M_PI;
  while (x_diff(3) < -M_PI) x_diff(3)+=2*M_PI;
  P_=P_ +weights_(i)*x_diff*x_diff.transpose();

}
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
int n_z=2;
  MatrixXd Zsig=MatrixXd(n_z,2*n_aug_+1);
  VectorXd z_pred=VectorXd(n_z);
   //Measurment of covraince matrix
   MatrixXd S=MatrixXd(n_z,n_z);

   //Transform sigma point into measurment space 
   for (int i=0;i<2*n_aug_+1;i++)
   {

     Zsig(0,i)=Xsig_pred_(0,i);
     Zsig(1,i)=Xsig_pred_(1,i);
   }
  //Mean prediticion measurment 
  z_pred(0.0);
  for (int i=0;i<2*n_aug_+1;i++)
   {
     z_pred=z_pred+weights_(i)*Zsig.col(i);
   }

   S.fill(0.0);
     for (int i=0;i<2*n_aug_+1;i++)
   {
  //Residual 
  VectorXd z_diff=Zsig.col(i)-z_pred;

  S=S+weights_(i)*z_diff*z_diff.transpose();
    }
//Add measurment noise covariance matrix
MatrixXd R=MatrixXd(n_z,n_z);
R<<std_laspx_*std_laspx_,0,
0,std_laspy_*std_laspy_;
S=S+R;

//Update 
MatrixXd Tc=MatrixXd(n_x_,n_z);
//Calcualte cross corelation matrix
Tc.fill(0.0);
  for (int i=0;i<2*n_aug_+1;i++)
   {
  //Residual 
  VectorXd z_diff=Zsig.col(i)-z_pred;
 
  //State difference 
   VectorXd x_diff=Xsig_pred_.col(i)-x_;
     //Angle normalization
  while (x_diff(3) > M_PI) x_diff(3)-=2*M_PI;
  while (x_diff(3) < -M_PI) x_diff(3)+=2*M_PI;
  Tc=Tc+weights_(i)*x_diff*z_diff.transpose();

    }
//Kalman gain 
MatrixXd K=Tc*S.inverse();

//Radar Measurment 
VectorXd z=meas_package.raw_measurements_;

//Residual 
  VectorXd z_diff=z-z_pred;
  //Angle normalization
  //while (z_diff(1) > M_PI) z_diff(1)-=2*M_PI;
  //while (z_diff(1) < -M_PI) z_diff(1)+=2*M_PI;

  //Update sate mean and covariance matrix
  x_=x_+K*z_diff ;
  P_ =P_ -K*S*K.transpose();

  //Calcualte NIS 
  NIS_Laser =z_diff.transpose() *S.inverse()*z_diff;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {

  int n_z=3;
  MatrixXd Zsig=MatrixXd(n_z,2*n_aug_+1);
  VectorXd z_pred=VectorXd(n_z);
   //Measurment of covraince matrix
   MatrixXd S=MatrixXd(n_z,n_z);

   //Transform sigma point into measurment space 
   Zsig.fill(0.0);
   for (int i=0;i<2*n_aug_+1;i++)
   {

     double p_x=Xsig_pred_(0,i);
     double p_y=Xsig_pred_(1,i);
     double v=Xsig_pred_(2,i);
     double yaw=Xsig_pred_(3,i);

     double v1=cos(yaw)*v;
     double v2=sin(yaw)*v;
     
     Zsig(0,i)=sqrt(p_x*p_x+p_y*p_y); //r
     Zsig(1,i)=atan2(p_y,p_x); //phi
      Zsig(2,i)=(p_x*v1+p_y*v2)/sqrt(p_x*p_x+p_y*p_y); //r_dot
   }
  //Mean prediticion measurment 
  z_pred.fill(0.0);
  for (int i=0;i<2*n_aug_+1;i++)
   {
     z_pred=z_pred+weights_(i)*Zsig.col(i);
   }

   S.fill(0.0);
     for (int i=0;i<2*n_aug_+1;i++)
   {
  //Residual 
  VectorXd z_diff=Zsig.col(i)-z_pred;
  //Angle normalization
  while (z_diff(1) > M_PI) z_diff(1)-=2*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1)+=2*M_PI;
  S=S+weights_(i)*z_diff*z_diff.transpose();
    }
//Add measurment noise covariance matrix
MatrixXd R=MatrixXd(n_z,n_z);
R<<std_radr_*std_radr_ ,0 ,0 ,
0,std_radphi_*std_radphi_,0,
0,0,std_radrd_*std_radrd_;
S=S+R;

//Update 
MatrixXd Tc=MatrixXd(n_x_,n_z);
//Calcualte cross corelation matrix
Tc.fill(0.0);
  for (int i=0;i<2*n_aug_+1;i++)
   {
  //Residual 
  VectorXd z_diff=Zsig.col(i)-z_pred;
  //Angle normalization
  while (z_diff(1) > M_PI) z_diff(1)-=2*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1)+=2*M_PI;

  //State difference 
   VectorXd x_diff=Xsig_pred_.col(i)-x_;
     //Angle normalization
  while (x_diff(3) > M_PI) x_diff(3)-=2*M_PI;
  while (x_diff(3) < -M_PI) x_diff(3)+=2*M_PI;
  Tc=Tc+weights_(i)*x_diff*z_diff.transpose();

    }
//Kalman gain 
MatrixXd K=Tc*S.inverse();

//Radar Measurment 
VectorXd z=meas_package.raw_measurements_;

//Residual 
  VectorXd z_diff=z-z_pred;
  //Angle normalization
  while (z_diff(1) > M_PI) z_diff(1)-=2*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1)+=2*M_PI;

  //Update sate mean and covariance matrix
  x_=x_+K*z_diff ;
  P_ =P_ -K*S*K.transpose();

  //Calcualte NIS 
  NIS_radar =z_diff.transpose() *S.inverse()*z_diff;

}