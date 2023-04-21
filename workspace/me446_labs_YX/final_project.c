#include "math.h"
#include "F28335Serial.h"

#define PI          3.1415926535897932384626433832795
#define TWOPI       6.283185307179586476925286766559
#define HALFPI      1.5707963267948966192313216916398
#define GRAV        9.81
#define NUMPOINTS 9

// These two offsets are only used in the main file user_CRSRobot.c  You just need to create them here and find the correct offset and then these offset will adjust the encoder readings
float offset_Enc2_rad = -22.95*PI/180.0;
float offset_Enc3_rad = 12.85*PI/180.0;

void cubic(float t, float *th1, float *th1d, float *th1dd);

typedef struct point_tag {
    float xb;
    float yb;
    float zb;
    float thetaz;
    int mode; //mode 0, stiff in all directions, mode 1 stiff only in z
} point;

point waypoints[NUMPOINTS] = {
                              {.15, 0, .43, 0, 0},
                              {.25, 0, .51, 0, 0}, //first point
                              {.04, .37, .25, 0, 0}, //slightly above hole
                              {.03, .35, .16, 0, 1}, //inside hole
                              {.03, .35, .4, 0, 1}, //raise outside hole
                              {.19, .12, .36, 0, 0}, //avoid obstacle
                              {.29, .10, .20, 0, 0}, //get to zigzag entrance
                              {.32, .05, .20, 0, 2}, //first zigzag point loose in x & y
                              {.23, .04, .20, 0, 2}, //second zigzag point
                              {.28, -.03, .19, 0, 2} //exit zigzag


};

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



// Your global varialbes.

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

float printtheta1motor = 0;
float printtheta2motor = 0;
float printtheta3motor = 0;

float x = 0;
float y = 0;
float z = 0;
float pi = PI;
float L1 = 0.254;
float L2 = 0.254;
float L3 = 0.254;

float invk_th1_dh = 0;
float invk_th2_dh = 0;
float invk_th3_dh = 0;
float invk_th1_m = 0;
float invk_th2_m = 0;
float invk_th3_m = 0;

float tau1print = 0;
float tau2print = 0;
float tau3print = 0;

//Friction compensation:
float errorBound = .05;

float minV1 = 0.1;
float steep1 = 3.6;
float viscousP1 = 0.1;
float coulombP1 = 0.3637;
float viscousN1 = 0.1;
float coulombN1 = -0.2948;

float minV2 = 0.05;
float steep2 = 9.;
float viscousP2 = 0.2;
float coulombP2 = 0.2;
float viscousN2 = 0.2;
float coulombN2 = -0.4;

float minV3 = 0.05;
float steep3 = 2;
float viscousP3 = 0.05;
float coulombP3 = 0.3;
float viscousN3 = 0.05;
float coulombN3 = -0.3;

float ffactor1 = 1;
float ffactor2 = 0.6;
float ffactor3 = 1;

float Theta1_old = 0;
float Omega1_old1 = 0;
float Omega1_old2 = 0;
float Omega1 = 0;

float Theta2_old = 0;
float Omega2_old1 = 0;
float Omega2_old2 = 0;
float Omega2 = 0;

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

float thetaz = 0; // atan(0.4/0.4) = pi/4 = 0.785 is the angle to rotate about z, make PDx gains weak
float thetax = 0;
float thetay = 0;

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

float x_a = 0;
float x_a_old = 0;
float y_a = 0;
float y_a_old = 0;
float z_a = 0;
float z_a_old = 0;

float x_ad = 0;
float x_ad_old1 = 0;
float x_ad_old2 = 0;

float y_ad = 0;
float y_ad_old1 = 0;
float y_ad_old2 = 0;

float z_ad = 0;
float z_ad_old1 = 0;
float z_ad_old2 = 0;

float x_des = 0.25;
float y_des = 0.25;
float z_des = 0.4;

float x_ddes = 0;
float y_ddes = 0;
float z_ddes = 0;

float Fx_des = 0;
float Fy_des = 0;
float Fz_des = 0;

float Fx_n = 0;
float Fy_n = 0;
float Fz_n = 0;
float Fx_w = 0;
float Fy_w = 0;
float Fz_w = 0;

float K_t = 6.0;


long total_time = 3000;
float progress = 0;

float x_des1 = 0;
float y_des1 = 0.4;
float z_des1 = 0.3;
float x_des2 = 0.4;
float y_des2 = 0.0;
float z_des2 = 0.3;

int targetWayPoint = 0;

// This function is called every 1 ms
void lab(float theta1motor,float theta2motor,float theta3motor,float *tau1,float *tau2,float *tau3, int error) {

//    //Forward Kinematics
    x = cos(theta1motor)*(L2*cos(pi/2 - theta2motor) + L3*cos(theta3motor));
    y = sin(theta1motor)*(L2*cos(pi/2 - theta2motor) + L3*cos(theta3motor));
    z = L1 + L2*sin(pi/2 - theta2motor) - L3*sin(theta3motor);

    invk_th1_dh = atan2(y, x);
    invk_th3_dh = acos(((z-L1)*(z-L1) + x*x + y*y - 2*L1*L1) / (2*L1*L1));
    invk_th2_dh = -atan2(z-L1, sqrt(x*x + y*y)) - invk_th3_dh/2;
//
    invk_th1_m = invk_th1_dh;
    invk_th2_m = invk_th2_dh + pi/2;
    invk_th3_m = invk_th3_dh + invk_th2_m - pi/2;


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

    //actual robot position and velocities
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



    //straight line
    progress = (float)(mycount % (total_time)) / (float)total_time;

    //rboot starting point
    x_des1 = waypoints[targetWayPoint].xb;
    y_des1 = waypoints[targetWayPoint].yb;
    z_des1 = waypoints[targetWayPoint].zb;

    x_des2 = waypoints[targetWayPoint + 1].xb;
    y_des2 = waypoints[targetWayPoint + 1].yb;
    z_des2 = waypoints[targetWayPoint + 1].zb;

    if(progress == 0){
        //progress = 0 ;
        targetWayPoint += 1;
    }

    if (targetWayPoint > NUMPOINTS) {
        //targetWayPoint = NUMPOINTS-1;
        //once all waypoints are reached
        x_des = .254;
        y_des = 0;
        z_des = .508;
    } else {
        x_des = x_des1 + progress * (x_des2 - x_des1);
        y_des = y_des1 + progress * (y_des2 - y_des1);
        z_des = z_des1 + progress * (z_des2 - z_des1);
    }




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

    float u_fric1;
    float u_fric2;
    float u_fric3;
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
   //tau1print = *tau1;
   //tau2print = *tau2;
   //tau3print = *tau3;
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


    float th1_m = 0;
    float th2_m = 0;
    float th3_m = 0;
    float th1_dh = th1_m;
    float th2_dh = th2_m - pi/2;
    float th3_dh = -th2_m + th3_m + pi/2;

    th1_m = invk_th1_m;
    th2_m = invk_th2_m;
    th3_m = invk_th3_m;

    th1_m = fmod(th1_m + pi, 2*pi) - pi;
    th2_m = fmod(th2_m + pi, 2*pi) - pi;
    th3_m = fmod(th3_m + pi, 2*pi) - pi;

    serial_printf(&SerialA, "%.2f (%.2f %.2f, %.2f) (%.2f %.2f, %.2f) (%.2f %.2f, %.2f)  \n\r", targetWayPoint, x,y,z, x_des1, y_des1, z_des1, x_des2, y_des2, z_des2);
    //serial_printf(&SerialA, "%.2f (%.2f %.2f, %.2f)    (%.2f %.2f, %.2f)  (%.2f, %.2f, %.2f)\n\r", t, printtheta1motor,printtheta2motor,printtheta3motor, th1_des,th2_des,th3_des,  tau1print,tau2print,tau3print);
}

