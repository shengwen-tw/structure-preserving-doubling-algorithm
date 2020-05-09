x1 = deg2rad(0);  %angle
x2 = 0;           %angular velocity
u = 0;            %torque, control input
g = 9.8;          %gravitational acceleration
l = 1;            %length of the rod [m]
k = 1;            %air drag coefficient
m = 0.5;          %mass [kg]

r = 2.4;          %gamma, parameter of Calyley transform
                  %SDA's author suggested the value between 2.1~2.6

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


%simulation time: sda_iteration_times * dt = ~10s
sda_iteration_times = 10000;
dt = 0.001;

%array for plotting
time_arr = zeros(1, sda_iteration_times);
x1_arr = zeros(1, sda_iteration_times);
x2_arr = zeros(1, sda_iteration_times);

tic();
for i = 1: sda_iteration_times
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

    G = B*inv(R)*transpose(B);
    H = transpose(C)*Q*C;
    A_trans = transpose(A);
    A_r = A - (r*I);
    
    %solve CARE with SDA
    sda_iteration_times = 0;
    I = eye(2);
    A_hat_last = I + 2*r*inv(A_r + G*inv(transpose(A_r))*H);
    G_hat_last = 2*r*inv(A_r)*G*inv(transpose(A_r) + H*inv(A_r)*G);
    H_hat_last = 2*r*inv(transpose(A_r) + H*inv(A_r)*G)*H*inv(A_r);

    while 1
        sda_iteration_times = sda_iteration_times + 1;
    
        %reduce redundent calculation by pre-calculating repeated terms
        inv_I_plus_H_G = inv(I + (H_hat_last * G_hat_last));
        transpose_A_hat_last = transpose(A_hat_last);
    
        %update
        A_hat_new = A_hat_last * inv(I + G_hat_last * H_hat_last) * A_hat_last;
        G_hat_new = G_hat_last + (A_hat_last * G_hat_last * inv_I_plus_H_G * transpose_A_hat_last);
        H_Hat_new = H_hat_last + (transpose_A_hat_last * inv_I_plus_H_G * H_hat_last * A_hat_last);
    
        %matrix norms
        norm_H_last = norm(H_hat_last);
        norm_H_now = norm(H_Hat_new);
    
        %prepare next iteration
        A_hat_last = A_hat_new;
        G_hat_last = G_hat_new;
        H_hat_last = H_Hat_new;
    
        %stop iteration if converged
        if abs(norm_H_now - norm_H_last) < 0.01
            break;
        end
    
        %disp(abs(norm_H_now - norm_H_last));
    end
    X = H_Hat_new;
    
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
title('eR');
xlabel('time [s]');
ylabel('angle [deg]');
subplot (2, 1, 2);
plot(time_arr, x2_arr);
xlabel('time [s]');
ylabel('angular acceleration [deg/s]');