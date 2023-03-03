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

float p1[]  = {0.25, -0.07, 0.17};
float p2[]  = {0.24, -0.03, 0.17};
float p3[]  = {0.25, -0.06, 0.17};
float p4[]  = {0.34, -0.05, 0.17};
float p5[]  = {0.34, -0.08, 0.17};
float p6[]  = {0.34, -0.03, 0.17};
float p7[]  = {0.28, -0.02, 0.21};
float p8[]  = {0.24, 0.03, 0.17};
float p9[]  = {0.24, -0.01, 0.17};
float p10[] = {0.34, -0.01, 0.17};
float p11[] = {0.33, 0.03, 0.17};
float p12[] = {0.28, 0.04, 0.23};
float p13[] = {0.24, 0.05, 0.17};
float p14[] = {0.28, 0.07, 0.17};
float p15[] = {0.24, 0.08, 0.17};
float p16[] = {0.28, 0.07, 0.17};
float p17[] = {0.33, 0.08, 0.17};
float p18[] = {0.30, 0.01, 0.23};
float p19[] = {0.25, -0.07, 0.17};
float p20[] = {0.25, -0.07, 0.17};
float p21[] = {0.25, -0.07, 0.17};

float* points[] = {p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, p21};

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

int pointID = 0;

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
float Kp2 = 50;
float Kd2 = 2.2;
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
float Kp3 = 50;
float Kd3 = 2;
float tau3print = 0;
float Ik3 = 0;
float ki3 = 400;

float errorBound = .05;
int period = 500;

// Assign these float to the values you would like to plot in Simulink
float Simulink_PlotVar1 = 0;
float Simulink_PlotVar2 = 0;
float Simulink_PlotVar3 = 0;
float Simulink_PlotVar4 = 0;


// This function is called every 1 ms
void lab(float theta1motor,float theta2motor,float theta3motor,float *tau1,float *tau2,float *tau3, int error) {

    //Forward Kinematics
    x = cos(theta1motor)*(L2*cos(pi/2 - theta2motor) + L3*cos(theta3motor));
    y = sin(theta1motor)*(L2*cos(pi/2 - theta2motor) + L3*cos(theta3motor));
    z = L1 + L2*sin(pi/2 - theta2motor) - L3*sin(theta3motor);




    //calculate cubic trajectory

    pointID = (int) ((mycount%(period*20)) / period); //increment pointID from 0 to 19 every 0.5s
    if (pointID > 19) {
        pointID = 0;
    }

    float prog = (mycount%period) / (float)(period);
    x = (1-prog)*points[pointID][0] + (prog)*points[pointID+1][0] + 0.1;
    y = (1-prog)*points[pointID][1] + (prog)*points[pointID+1][1];
    z = (1-prog)*points[pointID][2] + (prog)*points[pointID+1][2];

    invk_th1_dh = atan2(y, x);
    invk_th3_dh = acos(((z-L1)*(z-L1) + x*x + y*y - 2*L1*L1) / (2*L1*L1));
    invk_th2_dh = -atan2(z-L1, sqrt(x*x + y*y)) - invk_th3_dh/2;

    invk_th1_m = invk_th1_dh;
    invk_th2_m = invk_th2_dh + pi/2;
    invk_th3_m = invk_th3_dh + invk_th2_m - pi/2;

    th1_des = invk_th1_m;
    th2_des = invk_th2_m;
    th3_des = invk_th3_m;

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


    //PID
    *tau1 = Kp1*error1 + Kd1*error1d;
    *tau2 = Kp2*error2 + Kd2*error2d;
    *tau3 = Kp3*error3 + Kd3*error3d;


    //saturating torques
    if(fabs(*tau1) > 5 ||  fabs(error1) > errorBound){
        Ik1 = 0;
    } else {
        Ik1 += (error1 + error1old)/2.0 * 0.001;
    }
    if(fabs(*tau2) > 5 ||  fabs(error2) > errorBound){
        Ik2 = 0;
    }  else {
        Ik2 += (error2 + error2old)/2.0 * 0.001;
    }
    if(fabs(*tau3) > 5 ||   fabs(error3) > errorBound){
        Ik3 = 0;
    }  else {
        Ik3 += (error3 + error3old)/2.0 * 0.001;
    }

    *tau1 += ki1*Ik1;
    *tau2 += ki2*Ik2;
    *tau3 += ki3*Ik3;
//    *tau1 = 0;
//    *tau2 = 0;
//    *tau3 = 0;

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
   Simulink_PlotVar1 = theta1motor;
   Simulink_PlotVar2 = theta2motor;
   Simulink_PlotVar3 = theta3motor;
   Simulink_PlotVar4 = th2_des;

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

//    serial_printf(&SerialA, "%.2f (%.2f %.2f, %.2f)    (%.2f %.2f, %.2f)   (%.2f, %.2f, %.2f)  (%.2f, %.2f, %.2f)\n\r", t, printtheta1motor,printtheta2motor,printtheta3motor, th1_des,th2_des,th3_des, error1,error2,error3,  tau1print,tau2print,tau3print);
    serial_printf(&SerialA, "%d (%.2f %.2f, %.2f) (%.2f %.2f, %.2f) (%.2f %.2f, %.2f)  \n\r", pointID, x, y, z, th1_des,th2_des,th3_des, printtheta1motor,printtheta2motor,printtheta3motor);
}

