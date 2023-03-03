theta1_0 = 0;
theta1_1 = .5;
theta1_2 = 0;
theta1d_0 = 0;
theta1d_1 = 0;
theta1d_2 = 0;

t0 = 0;
t1 = 1;
t2 = 2;

A01 = [[1 t0 t0^2 t0^3];
      [0 1 2*t0 3*t0^2];
      [1 t1 t1^2 t1^3];
      [0 1 2*t1 3*t1^2];];
b01 = [theta1_0; theta1d_0; theta1_1; theta1d_1];
x01 = A01\b01


A12 = [[1 t1 t1^2 t1^3];
      [0 1 2*t1 3*t1^2];
      [1 t2 t2^2 t2^3];
      [0 1 2*t2 3*t2^2];];
b12 = [theta1_1; theta1d_1; theta1_2; theta1d_2];
x12 = A12\b12


%state = [theta1_0; theta1d_0; theta1dd_0]


t = linspace(0,3);
states = Cubic(x01, x12, 0)


function state = Cubic(x01, x12, t)
if(t < 1)
    a0=x01(1);
    a1=x01(2);
    a2=x01(3);
    a3=x01(4);
else
    a0=x12(1);
    a1=x12(2);
    a2=x12(3);
    a3=x12(4);
end
th1 = a0 + a1*t + a2*t^2 + a3*t^3;
th1d = a1 + 2*a2*t + 3*a3*t^2;
th1dd = 2*a2 + 6*a3;
if(t > 0)
    state = [0 0 0];
else
state = [th1 th1d th1dd];
end
end


