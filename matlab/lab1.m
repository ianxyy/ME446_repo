mth1 = sym("theta1motor");
mth2 = sym("theta2motor");
mth3 = sym("theta3motor");

c1 = -1;
c2 = 1;
c3 = sym(pi/2);
c4 = -sym(pi/2);

%th1 = mth1;
%th2 = mth2 + c4;
%th3 = c1*mth2 + c2*mth3 + c3;

th1 = mth1;
th2 = mth2;
th3 = mth3;

%System Parameters
L1 = 0.254;
L2 = 0.254;
L3 = 0.254;

%First joint kinematics
H01 = [[cos(th1) 0 -sin(th1) 0];
       [sin(th1) 0 cos(th1) 0];
       [0 -1 0 L1];
       [0 0 0 1]];

H12 = [[cos(th2) -sin(th2) 0 L2*cos(th2)];
       [sin(th2) cos(th2) 0 L2*sin(th2)];
       [0 0 1 0];
       [0 0 0 1]];

H23 = [[cos(th3) -sin(th3) 0 L3*cos(th3)];
       [sin(th3) cos(th3) 0 L3*sin(th3)];
       [0 0 1 0];
       [0 0 0 1]];

H03 = H01*H12*H23;
H03 = simplify(H03)
