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

float J1 = 0.0167;
float J2 = 0.03;
float J3 = 0.0128;

float invk_th1_dh = 0;
float invk_th2_dh = 0;
float invk_th3_dh = 0;
float invk_th1_m = 0;
float invk_th2_m = 0;
float invk_th3_m = 0;

float th1_des = 0;
float th1d_des = 0;
float th1dd_des = 0;
float error1 = 0;
float error1d = 0;
float error1old = 0;
float Theta1_old = 0;
float Omega1_old1 = 0;
float Omega1_old2 = 0;
float Omega1 = 0;
float Kp1 = 50;
float Kd1 = 2;
float tau1print = 0;
float Ik1 = 0;
float ki1 = 400;

float t = 0;

float th2_des = 0;
float th2d_des = 0;
float th2dd_des = 0;
float error2 = 0;
float error2d = 0;
float error2old = 0;
float Theta2_old = 0;
float Omega2_old1 = 0;
float Omega2_old2 = 0;
float Omega2 = 0;
float Kp2 = 2000;
float Kd2 = 100;
float tau2print = 0;
float Ik2 = 0;
float ki2 = 450;


float th3_des = 0;
float th3d_des = 0;
float th3dd_des = 0;
float error3 = 0;
float error3d = 0;
float error3old = 0;
float Theta3_old = 0;
float Omega3_old1 = 0;
float Omega3_old2 = 0;
float Omega3 = 0;
float Kp3 = 3000;
float Kd3 = 60;
float tau3print = 0;
float Ik3 = 0;
float ki3 = 400;

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

//used to calculate DCG matrices
float p1 = 0.0300;
float p2 = 0.0128;
float p3 = 0.0076;
float p4 = 0.0753;
float p5 = 0.0298;


// Assign these float to the values you would like to plot in Simulink
float Simulink_PlotVar1 = 0;
float Simulink_PlotVar2 = 0;
float Simulink_PlotVar3 = 0;
float Simulink_PlotVar4 = 0;
float Simulink_PlotVar5 = 0;

int doCubic = 0;

float traj2_low = 0;
float traj2_high = PI/6;
float traj3_low = 0;
float traj3_high = PI/6;
float th2_step = .25;
float th3_step = .3;

int doInverseDynamics = 0;
//PD gains without inverse dynamics:
float Kp2_PD = 60;
float Kd2_PD = 5;
float Kp3_PD = 60;
float Kd3_PD = 5;


// This function is called every 1 ms
void lab(float theta1motor,float theta2motor,float theta3motor,float *tau1,float *tau2,float *tau3, int error) {

    //Forward Kinematics
    x = cos(theta1motor)*(L2*cos(pi/2 - theta2motor) + L3*cos(theta3motor));
    y = sin(theta1motor)*(L2*cos(pi/2 - theta2motor) + L3*cos(theta3motor));
    z = L1 + L2*sin(pi/2 - theta2motor) - L3*sin(theta3motor);

    invk_th1_dh = atan2(y, x);
    invk_th3_dh = acos(((z-L1)*(z-L1) + x*x + y*y - 2*L1*L1) / (2*L1*L1));
    invk_th2_dh = -atan2(z-L1, sqrt(x*x + y*y)) - invk_th3_dh/2;

    invk_th1_m = invk_th1_dh;
    invk_th2_m = invk_th2_dh + pi/2;
    invk_th3_m = invk_th3_dh + invk_th2_m - pi/2;


    //calculate cubic trajectory
    if(doCubic){
        t = (mycount%3000) / 1000.0;
        float th_cubic = 0;
        float thd_cubic = 0;
        float thdd_cubic = 0;
        cubic(t, &th_cubic, &thd_cubic, &thdd_cubic);
        th1_des = th_cubic;
        th2_des = th_cubic;
        th3_des = th_cubic;
        th1d_des = thd_cubic;
        th2d_des = thd_cubic;
        th3d_des = thd_cubic;
        th1dd_des = thdd_cubic;
        th2dd_des = thdd_cubic;
        th3dd_des = thdd_cubic;
    }else{
        //square wave between 0 and pi/6
        if(mycount % 3000 == 0){
            if(th2_step  > .4){
                th2_step = .25;
            }else{
                th2_step = .85;
            }
            if(th3_step < 0){
                th3_step = .3;
            }else{
                th3_step = -.3;
            }
        }
        implement_discrete_tf(&trajectory2, th2_step, &th2_des, &th2d_des, &th2dd_des);
        th1_des = th2_des;
        th1d_des = th2d_des;
        th1dd_des = th2dd_des;
        implement_discrete_tf(&trajectory3, th3_step, &th3_des, &th3d_des, &th3dd_des);
    }

    //calc states

    //get th2dot from IIR filter
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

    error1 = th1_des - theta1motor;
    error1d = th1d_des - Omega1;
    error2 = th2_des - theta2motor;
    error2d = th2d_des - Omega2;
    error3 = th3_des - theta3motor;
    error3d = th3d_des - Omega3;

    //calculate desired accelerations a
    float ath2 = th2dd_des + Kp2*(error2) + Kd2*(error2d);
    float ath3 = th3dd_des + Kp3*(error3) + Kd3*(error3d);

    //calculate DCG matrices
    float sin32 = sin(theta3motor - theta2motor);
    float cos32 = cos(theta3motor - theta2motor);

    float innerD2 = p1*ath2 - p3*sin32*ath3;
    float innerD3 = -p3*ath2*sin32 + p2*ath3;
    float innerC2 = -Omega3*p3*cos32*Omega3;
    float innerC3 = Omega2*p3*cos32*Omega2;
    float innerG2 = -p4*GRAV*sin(theta2motor);
    float innerG3 = -p5*GRAV*cos(theta3motor);

    //PID
    *tau1 = J1*th1dd_des + Kp1*error1 + Kd1*error1d;

    if(doInverseDynamics){
        *tau2 = innerD2 + innerC2 + innerG2;
        *tau3 = innerD3 + innerC3 + innerG3;
    }else{
        //do PD + feedforward
        *tau2 = J2*th2dd_des + Kp2_PD*error2 + Kd2_PD*error2d;
        *tau3 = J3*th3dd_des + Kp3_PD*error3 + Kd3_PD*error3d;
    }



    //saturating torques
//    if(fabs(*tau1) > 5 ||  fabs(error1) > errorBound){
//        Ik1 = 0;
//    } else {
//        Ik1 += (error1 + error1old)/2.0 * 0.001;
//    }
//    if(fabs(*tau2) > 5 ||  fabs(error2) > errorBound){
//        Ik2 = 0;
//    }  else {
//        Ik2 += (error2 + error2old)/2.0 * 0.001;
//    }
//    if(fabs(*tau3) > 5 ||   fabs(error3) > errorBound){
//        Ik3 = 0;
//    }  else {
//        Ik3 += (error3 + error3old)/2.0 * 0.001;
//    }

//    *tau1 += ki1*Ik1;
//    *tau2 += ki2*Ik2;
//    *tau3 += ki3*Ik3;

    //friction compensation
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
    *tau1 += u_fric1;
    *tau2 += u_fric2;
    *tau3 += u_fric3;

    error1old = error1;
    error2old = error2;
    error3old = error3;



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
   Simulink_PlotVar1 = theta2motor;
   Simulink_PlotVar2 = th2_des;
   Simulink_PlotVar3 = theta3motor;
   Simulink_PlotVar4 = th3_des;

   //LED
   if ((mycount%500)==0) {
       UARTprint = 1;
       GpioDataRegs.GPBTOGGLE.bit.GPIO34 = 1; // Blink LED on Control Card
       GpioDataRegs.GPBTOGGLE.bit.GPIO60 = 1; // Blink LED on Emergency Stop Box
   }


    mycount++;
}

void cubic(float t, float *th1, float *th1d, float *th1dd){
   float a0 = 0;
   float a1 = 0;
   float a2 = 0;
   float a3 = 0;
    if(t < 1){
       a0=0;
       a1=0;
       a2=1.5;
       a3=-1;
   }else{
       a0=-2;
       a1=6;
       a2=-4.5;
       a3=1;
   }
   *th1 = a0 + a1*t + a2*t*t + a3*t*t*t;
   *th1d = a1 + 2*a2*t + 3*a3*t*t;
   *th1dd = 2*a2 + 6*a3;
   if(t < 0 || t > 2){
       *th1 = 0;
       *th1d = 0;
       *th1dd = 0;
   }
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

    serial_printf(&SerialA, "%.2f (%.2f %.2f, %.2f)    (%.2f %.2f, %.2f)   (%.2f, %.2f, %.2f)  (%.2f, %.2f, %.2f)\n\r", t, printtheta1motor,printtheta2motor,printtheta3motor, th1_des,th2_des,th3_des, error1,error2,error3,  tau1print,tau2print,tau3print);
}

