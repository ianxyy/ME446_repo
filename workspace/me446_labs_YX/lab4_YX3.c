#include "math.h"
#include "F28335Serial.h"

#define PI          3.1415926535897932384626433832795
#define TWOPI       6.283185307179586476925286766559
#define HALFPI      1.5707963267948966192313216916398
#define GRAV        9.81

// These two offsets are only used in the main file user_CRSRobot.c  You just need to create them here and find the correct offset and then these offset will adjust the encoder readings
float offset_Enc2_rad = -22.95*PI/180.0;
float offset_Enc3_rad = 12.85*PI/180.0;

void cubic(float t, float *th1, float *th1d, float *th1dd);



typedef struct steptraj_s {
    long double b[4];
    long double a[4];
    long double xk[4];
    long double yk[4];
    float qd_old;
    float qddot_old;
    int size;
} steptraj_t;

steptraj_t trajectory2 = {9.7059014792764445e-07L,2.9117704437829331e-06L,2.9117704437829331e-06L,9.7059014792764445e-07L,
                        1.0000000000000000e+00L,-2.9405940594059405e+00L,2.8823644740711698e+00L,-9.4176264994404557e-01L,
                        0,0,0,0,
                        0,0,0,0,
                        0,
                        0,
                        4};

steptraj_t trajectory3 = {9.7059014792764445e-07L,2.9117704437829331e-06L,2.9117704437829331e-06L,9.7059014792764445e-07L,
                        1.0000000000000000e+00L,-2.9405940594059405e+00L,2.8823644740711698e+00L,-9.4176264994404557e-01L,
                        0,0,0,0,
                        0,0,0,0,
                        0,
                        0,
                        4};

// this function must be called every 1ms.
void implement_discrete_tf(steptraj_t *traj, float step, float *qd, float *qd_dot, float *qd_ddot) {
    int i = 0;

    traj->xk[0] = step;
    traj->yk[0] = traj->b[0]*traj->xk[0];
    for (i = 1;i<traj->size;i++) {
        traj->yk[0] = traj->yk[0] + traj->b[i]*traj->xk[i] - traj->a[i]*traj->yk[i];
    }

    for (i = (traj->size-1);i>0;i--) {
        traj->xk[i] = traj->xk[i-1];
        traj->yk[i] = traj->yk[i-1];
    }

    *qd = traj->yk[0];
    *qd_dot = (*qd - traj->qd_old)*1000;  //0.001 sample period
    *qd_ddot = (*qd_dot - traj->qddot_old)*1000;

    traj->qd_old = *qd;
    traj->qddot_old = *qd_dot;
}



// global variables.

//counter used in lab() function for determining when to do certain tasks
long mycount = 0;

#pragma DATA_SECTION(whattoprint, ".my_vars")
float whattoprint = 0.0;

#pragma DATA_SECTION(theta1array, ".my_arrs")
float theta1array[100];

#pragma DATA_SECTION(print2, ".my_vars")
float print2 = 12345678.9;

#pragma DATA_SECTION(theta2array, ".my_arrs")
float theta2array[100];

long arrayindex = 0;
int UARTprint = 0;

//variables used for printing out our motor angles
float printtheta1motor = 0;
float printtheta2motor = 0;
float printtheta3motor = 0;

//variables to define robot link lengths
float pi = PI;
float L1 = 0.254;
float L2 = 0.254;
float L3 = 0.254;

//variables used to print out each motor's torque
float tau1print = 0;
float tau2print = 0;
float tau3print = 0;

//Friction compensation variables corresponding to each joint and velocity conditions

//Motor 1 friction compensation where we define the coefficients for positive and negative coulomb & viscous friction coefficients
// as well as the minimum velocity required to switch between coulomb and viscous friciton. steep1 refers to steep slope region near 0 velocity
float minV1 = 0.1;
float steep1 = 3.6;
float viscousP1 = 0.1;
float coulombP1 = 0.3637;
float viscousN1 = 0.1;
float coulombN1 = -0.2948;

//Motor 2 friction compensation where the variables have the same definition as those for motor 1
float minV2 = 0.05;
float steep2 = 9.;
float viscousP2 = 0.2;
float coulombP2 = 0.2;
float viscousN2 = 0.2;
float coulombN2 = -0.4;

//motor 3 friciton compensation
float minV3 = 0.05;
float steep3 = 2;
float viscousP3 = 0.05;
float coulombP3 = 0.3;
float viscousN3 = 0.05;
float coulombN3 = -0.3;

//variable to scale the amount of friction compensation we want affecting our joints, 1 is full friction compensation & .6 is 60% friction compensation
float ffactor1 = 1;
float ffactor2 = 0.6;
float ffactor3 = 1;

//variables used to calculate actual filtered velocity (Omega) with IIR filter for each joint

//Omega is actual velocity we calculate, Theta is our actual motor angle
float Theta1_old = 0;
float Omega1_old1 = 0;
float Omega1_old2 = 0;
float Omega1 = 0;

//variables used to calculate filtered velocity for joint 2
float Theta2_old = 0;
float Omega2_old1 = 0;
float Omega2_old2 = 0;
float Omega2 = 0;

//variables used to calculate filtered velocity for joint 3
float Theta3_old = 0;
float Omega3_old1 = 0;
float Omega3_old2 = 0;
float Omega3 = 0;



// Assign these float to the values you would like to plot in Simulink
float Simulink_PlotVar1 = 0;
float Simulink_PlotVar2 = 0;
float Simulink_PlotVar3 = 0;
float Simulink_PlotVar4 = 0;



//PD gains without inverse dynamics:
float Kp_x = 500;
float Kd_x = 25;
float Kp_y = 500;
float Kd_y = 25;
float Kp_z = 500;
float Kd_z = 25;

//variables used to define how our impedance control changes due to rotations along world axes
float thetaz = 0; // atan(0.4/0.4) = pi/4 = 0.785 is the angle to rotate about z, make PDx gains weak
float thetax = 0;
float thetay = 0;

//variables used to define task space control
float cosq1 = 0;
float sinq1 = 0;
float cosq2 = 0;
float sinq2 = 0;
float cosq3 = 0;
float sinq3 = 0;
float JT_11 = 0;
float JT_12 = 0;
float JT_13 = 0;
float JT_21 = 0;
float JT_22 = 0;
float JT_23 = 0;
float JT_31 = 0;
float JT_32 = 0;
float JT_33 = 0;
float cosz = 0;
float sinz = 0;
float cosx = 0;
float sinx = 0;
float cosy = 0;
float siny = 0;

//variables used to define entries to our rotation matricies for transforming impedance control axes
float R11 = 0;
float R12 = 0;
float R13 = 0;
float R21 = 0;
float R22 = 0;
float R23 = 0;
float R31 = 0;
float R32 = 0;
float R33 = 0;
float RT11 = 0;
float RT12 = 0;
float RT13 = 0;
float RT21 = 0;
float RT22 = 0;
float RT23 = 0;
float RT31 = 0;
float RT32 = 0;
float RT33 = 0;

//variables used to calculate actual task space coordinates, x_a, y_a, z_a
float x_a = 0;
float x_a_old = 0;
float y_a = 0;
float y_a_old = 0;
float z_a = 0;
float z_a_old = 0;

//variables used to store actual task space velocities which are found via IIR filter
float x_ad = 0;
float x_ad_old1 = 0;
float x_ad_old2 = 0;

//variables used to find task space y velocity
float y_ad = 0;
float y_ad_old1 = 0;
float y_ad_old2 = 0;

//variables used to store task space z velocity
float z_ad = 0;
float z_ad_old1 = 0;
float z_ad_old2 = 0;

//variables used to define our desired task space position
float x_des = 0.25;
float y_des = 0.25;
float z_des = 0.4;

//variables used to define our desired task space velocities
float x_ddes = 0;
float y_ddes = 0;
float z_ddes = 0;

//variables used to define our desired task space forces
float Fx_des = 0;
float Fy_des = 0;
float Fz_des = 0;

//variables used to store forces in the world frame and their values in a transformed frame, n
float Fx_n = 0;
float Fy_n = 0;
float Fz_n = 0;
float Fx_w = 0;
float Fy_w = 0;
float Fz_w = 0;

//scale factor for our commanded forces
float K_t = 6.0;

//variables used to follow a line trajectory

//total time allowed for our trajectory to be completed in
long total_time = 2000;

//variable that acts as a 'percentage' to determine how much of the way we are to our desired position
float progress = 0;

//task space coordinates that define our two desired endpoints for following a line trajectory
float x_des1 = 0;
float y_des1 = 0.4;
float z_des1 = 0.3;
float x_des2 = 0.4;
float y_des2 = 0.0;
float z_des2 = 0.3;

// This function is called every 1 ms
void lab(float theta1motor,float theta2motor,float theta3motor,float *tau1,float *tau2,float *tau3, int error) {
    // Rotation zxy and its Transpose
    cosz = cos(thetaz);
    sinz = sin(thetaz);
    cosx = cos(thetax);
    sinx = sin(thetax);
    cosy = cos(thetay);
    siny = sin(thetay);
    //R changes from N frame to world frame
    //RT changes from world frame to N frame
    RT11 = R11 = cosz*cosy-sinz*sinx*siny;
    RT21 = R12 = -sinz*cosx;
    RT31 = R13 = cosz*siny+sinz*sinx*cosy;
    RT12 = R21 = sinz*cosy+cosz*sinx*siny;
    RT22 = R22 = cosz*cosx;
    RT32 = R23 = sinz*siny-cosz*sinx*cosy;
    RT13 = R31 = -cosx*siny;
    RT23 = R32 = sinx;
    RT33 = R33 = cosx*cosy;
    // Jacobian Transpose
    cosq1 = cos(theta1motor);
    sinq1 = sin(theta1motor);
    cosq2 = cos(theta2motor);
    sinq2 = sin(theta2motor);
    cosq3 = cos(theta3motor);
    sinq3 = sin(theta3motor);
    JT_11 = -0.254*sinq1*(cosq3 + sinq2);
    JT_12 = 0.254*cosq1*(cosq3 + sinq2);
    JT_13 = 0;
    JT_21 = 0.254*cosq1*(cosq2 - sinq3);
    JT_22 = 0.254*sinq1*(cosq2 - sinq3);
    JT_23 = -0.254*(cosq3 + sinq2);
    JT_31 = -0.254*cosq1*sinq3;
    JT_32 = -0.254*sinq1*sinq3;
    JT_33 = -0.254*cosq3;

    /**
    We first calculate our current states in the task space which are our actual position and actual velocities. For finding our actual velocities we take
    a measurement of our task space coordinates in two consequtive lab() iterations and approximate the velocity as the difference in position over the 
    elapsed time of 1ms, as that is how often the lab() function runs. We then smooth our velocities by implementing and IIR filter.
    */
    
    x_a = 0.254*cosq1*(cosq3+sinq2);
    y_a = 0.254*sinq1*(cosq3+sinq2);
    z_a = 0.254*(1+cosq2-sinq3);

    x_ad = (x_a - x_a_old)/0.001;
    x_ad = (x_ad + x_ad_old1 + x_ad_old2)/3.0;
    x_a_old = x_a;
    x_ad_old2 = x_ad_old1;
    x_ad_old1 = x_ad;

    y_ad = (y_a - y_a_old)/0.001;
    y_ad = (y_ad + y_ad_old1 + y_ad_old2)/3.0;
    y_a_old = y_a;
    y_ad_old2 = y_ad_old1;
    y_ad_old1 = y_ad;

    z_ad = (z_a - z_a_old)/0.001;
    z_ad = (z_ad + z_ad_old1 + z_ad_old2)/3.0;
    z_a_old = z_a;
    z_ad_old2 = z_ad_old1;
    z_ad_old1 = z_ad;

    /**
    After finding our current state, we are able to then follow our desired trajectory which is a line between two points. In order to do this we defined a
    progress variable which is a multiplier that takes values between 0 and 1. When the progress variable is at 0, then our desired position will correspond
    to our first desired point. As progress approaches the value of 1, then our desired position will become our second desired point. The function that 
    gives us our desired position based on our progress and the two endpoints is then:
    
        desired_position = endpoint_1 + progress * (endpoint2 - endpoint1)
    
    
    Instead of repeating this twice and flipping our desired points, we instead allow our progress variable to go above 1, but now subtract a value of 2 from
    it since we will be starting at our second desired point, therefore 2 - progress = 2 - 1 = 1, so we start at our second point and as progress continues to
    increase up to the value of 2, then we will have 2 - progess = 2 - 2 = 0, so we will then be at our first desired point.
    */
    
    //define progress variable to range from 0 to 2 and repeat indefinately to repeat our line trajectory back and forth
    progress = (float)(mycount % (2*total_time)) / (float)total_time;
    if(progress > 1){
        progress = 2 - progress;
    }

    //update our desired task space coordinates based on our current progress
    x_des = x_des1 + progress * (x_des2 - x_des1);
    y_des = y_des1 + progress * (y_des2 - y_des1);
    z_des = z_des1 + progress * (z_des2 - z_des1);


    /**
        After calculating our desired position based off our progress in the trajectory, we can choose to rotate our coordinate axes that define our impedance
        control in order to make the stiffness in any direction we want be either weaker or stronger. 
    */
    
    //converts world errors to N frame errors
    float x_error = RT11*(x_des - x_a) + RT12*(y_des - y_a) + RT13*(z_des - z_a);
    float y_error = RT21*(x_des - x_a) + RT22*(y_des - y_a) + RT23*(z_des - z_a);
    float z_error = RT31*(x_des - x_a) + RT32*(y_des - y_a) + RT33*(z_des - z_a);

    float x_derror = RT11*(x_ddes - x_ad) + RT12*(y_ddes - y_ad) + RT13*(z_ddes - z_ad);
    float y_derror = RT21*(x_ddes - x_ad) + RT22*(y_ddes - y_ad) + RT23*(z_ddes - z_ad);
    float z_derror = RT31*(x_ddes - x_ad) + RT32*(y_ddes - y_ad) + RT33*(z_ddes - z_ad);

    Fx_n = Kp_x * x_error + Kd_x * x_derror;
    Fy_n = Kp_y * y_error + Kd_y * y_derror;
    Fz_n = Kp_z * z_error + Kd_z * z_derror;

    Fx_w = R11 * Fx_n + R12 * Fy_n + R13 * Fz_n       + Fx_des/K_t;
    Fy_w = R21 * Fx_n + R22 * Fy_n + R23 * Fz_n       + Fy_des/K_t;
    Fz_w = R31 * Fx_n + R32 * Fy_n + R33 * Fz_n       + Fz_des/K_t;

    //task space PD
     *tau1 = JT_11 * Fx_w + JT_12 * Fy_w + JT_13 * Fz_w;
     *tau2 = JT_21 * Fx_w + JT_22 * Fy_w + JT_23 * Fz_w;
     *tau3 = JT_31 * Fx_w + JT_32 * Fy_w + JT_33 * Fz_w;

    /**
        After creating our initial task space PD controller we then wanted to add friction compensation, and in order to do so we first need to calculate
        our joint's velocities using an IIR filter. Since we know that our lab() function runs every 1ms we can approximate the velocity as the difference
        of two consecutive motor angle measurements divided by the elapsed time of 1ms. This would be very noisy however, so we save the approximated velocity
        and average it out with the previous two measurements in order to have a smoother velocity for each joint.
    */
    
    //friction compensation and calc thetadots
     Omega1 = (theta1motor - Theta1_old)/0.001;
     Omega1 = (Omega1 + Omega1_old1 + Omega1_old2)/3.0;
     Theta1_old = theta1motor;
     Omega1_old2 = Omega1_old1;
     Omega1_old1 = Omega1;

     Omega2 = (theta2motor - Theta2_old)/0.001;
     Omega2 = (Omega2 + Omega2_old1 + Omega2_old2)/3.0;
     Theta2_old = theta2motor;
     Omega2_old2 = Omega2_old1;
     Omega2_old1 = Omega2;

     Omega3 = (theta3motor - Theta3_old)/0.001;
     Omega3 = (Omega3 + Omega3_old1 + Omega3_old2)/3.0;
     Theta3_old = theta3motor;
     Omega3_old2 = Omega3_old1;
     Omega3_old1 = Omega3;
    
    /**
        After calculating our filtered joint velocities, Omega_i, we then go on to calculate our friction compensation for each joint. We calculate the friction 
        compensation for joint i by first checking what region on the friction vs velocity graph it lies on. There are 3 regions we have defined:
        
            1) Joint velocity is greater than a minimum velocity due to the friction being nonlinear near zero
            2) Joint velocity is less than a minimum velocity again due to the friction being nonlinear near zero
            3) Joint velocity is in the regious between the other regions typically near zero
        
        Once we have determined which regious we are operating in, we are able to create a friction compensating force for joint i in the form of:
            
            u_friction_i = viscous_friction_coefficient * joint_angular_velocity + coulomb_friction
       
       if the angular velocity is outside of the operating region near zero, and if the angular velocity, Omega, is near zero then we apply a force of:
            
            u_friction_i = coefficient * joint_angular_velocity
       
       Since each of our joints can be moving at a different angular velocity, we have to find the operating region of each motor and then define 
       the corresponding friction compensation forces, u_fric1, u_fric2, u_fric3, one by one.
       
       Once we have the friction compensation forces we multiply them by a scale factor incase we want to have less or more friction compensation for each
       joint, and lastly we add this to our current control law.
    **/

    float u_fric1;
    float u_fric2;
    float u_fric3;
    
    //define friction compensation depending on operating region for each joint
    if (Omega1 > minV1) {
            u_fric1 = viscousP1*Omega1 + coulombP1;
        } else if (Omega1 < -minV1) {
            u_fric1 = viscousN1*Omega1 + coulombN1;
        } else {
            u_fric1 = steep1*Omega1;
    }
    
    if (Omega2 > minV2) {
            u_fric2 = viscousP2*Omega2 + coulombP2;
        } else if (Omega2 < -minV2) {
            u_fric2 = viscousN2*Omega2 + coulombN2;
        } else {
            u_fric2 = steep2*Omega2;
        }
    
    if (Omega3 > minV3) {
            u_fric3 = viscousP3*Omega3 + coulombP3;
        } else if (Omega3 < -minV3) {
            u_fric3 = viscousN3*Omega3 + coulombN3;
        } else {
            u_fric3 = steep3*Omega3;
    }

    //add scaled friction to current control effort
    *tau1 += u_fric1 * ffactor1;
    *tau2 += u_fric2 * ffactor2;
    *tau3 += u_fric3 * ffactor3;


    //save motor angles for plotting
   if ((mycount%50)==0) {
       theta1array[arrayindex] = theta1motor;
       theta2array[arrayindex] = theta2motor;
       if (arrayindex >= 100) {
           arrayindex = 0;
       } else {
           arrayindex++;
       }
   }

   printtheta1motor = theta1motor*180/PI;
   printtheta2motor = theta2motor*180/PI;
   printtheta3motor = theta3motor*180/PI;
   tau1print = *tau1;
   tau2print = *tau2;
   tau3print = *tau3;
   Simulink_PlotVar1 = x_a;
   Simulink_PlotVar2 = x_des;
   Simulink_PlotVar3 = y_a;
   Simulink_PlotVar4 = y_des;

   //LED
   if ((mycount%500)==0) {
       UARTprint = 1;
       GpioDataRegs.GPBTOGGLE.bit.GPIO34 = 1; // Blink LED on Control Card
       GpioDataRegs.GPBTOGGLE.bit.GPIO60 = 1; // Blink LED on Emergency Stop Box
   }


    mycount++;
}


void printing(void){


//    float th1_m = 0;
//    float th2_m = 0;
//    float th3_m = 0;
//    float th1_dh = th1_m;
//    float th2_dh = th2_m - pi/2;
//    float th3_dh = -th2_m + th3_m + pi/2;

//    th1_m = invk_th1_m;
//    th2_m = invk_th2_m;
//    th3_m = invk_th3_m;

    //th1_m = fmod(th1_m + pi, 2*pi) - pi;
    //th2_m = fmod(th2_m + pi, 2*pi) - pi;
    //th3_m = fmod(th3_m + pi, 2*pi) - pi;

//    serial_printf(&SerialA, "%.2f (%.2f %.2f, %.2f)    (%.2f %.2f, %.2f)   (%.2f, %.2f, %.2f)  (%.2f, %.2f, %.2f)\n\r", t, printtheta1motor,printtheta2motor,printtheta3motor, th1_des,th2_des,th3_des, error1,error2,error3,  tau1print,tau2print,tau3print);
}

