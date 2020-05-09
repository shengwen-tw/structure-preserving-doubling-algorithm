x1 = deg2rad(0);  %angle
x2 = 0;           %angular velocity
u = 0;            %torque, control input
g = 9.8;          %gravitational acceleration
l = 1;            %length of the rod [m]
k = 1;            %air drag coefficient
m = 0.5;          %mass [kg]

%parameter for full state feedback control
K_full_state_fb = [100 1]; %k1: angle control gain,
                           %k2: angular velocity control gain

%parameters for H2 control
%Q: state penalty
Q = [3000 0
     0 10];
%R: control penalty
R = [0.1];

%controller setpoint
x1_d = deg2rad(100); %desired angle
x2_d = deg2rad(0); %desired angular velocity
x_d = [x1_d; x2_d];

%simulation time: ITERATION_TIMES * dt = ~10s
ITERATION_TIMES = 10000;
dt = 0.001;

%array for plotting
time_arr = zeros(1, ITERATION_TIMES);
x1_arr = zeros(1, ITERATION_TIMES);
x2_arr = zeros(1, ITERATION_TIMES);

tic();
for i = 1: ITERATION_TIMES
    %update system dynamics
    x1_dot = x2;
    x2_dot = (-g/l)*sin(x1) - (k/m)*x2 + (1/(m*l*l))*u;
    x1 = x1 + x1_dot * dt;
    x2 = x2 + x2_dot * dt;
    x = [x1; x2];

    %linear quadratic regulator
    a11 = 0;
    a12 = 1;
    a21 = (-g/l)*cos(x1);
    a22 = -k/m;
    A = [a11 a12;
         a21 a22];
    
    b1 = 0;
    b2 = 1/(m*l*l);
    B = [b1;
         b2];
    
    %solve CARE
    C = eye(2); %measurement matrix is identical, every states are measurable
    H = transpose(C) * Q * C;
    [X, L, G] = care(A, B, H, R);
    
    %calculate optimal feedback control gain and control input
    K_lqr = inv(R) * transpose(B) * X;
    u = -K_lqr*(x-x_d);
    
    %full state feedback control
    %u = -K_full_state_fb*(x-x_d);
    
    %record datas for plotting
    x1_arr(i) = rad2deg(x1);
    x2_arr(i) = rad2deg(x2);
    time_arr(i) = i * dt;
end
toc();

figure(1);
subplot (2, 1, 1);
plot(time_arr, x1_arr);
title('Pendulum LQR Control');
xlabel('time [s]');
ylabel('angle [deg]');
subplot (2, 1, 2);
plot(time_arr, x2_arr);
xlabel('time [s]');
ylabel('angular velocity [deg/s]');