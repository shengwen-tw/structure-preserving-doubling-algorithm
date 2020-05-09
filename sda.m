%inputs
A = [-3 2;
      1 1];
B = [0; 1];
C = [1 -1];
R = [1];
Q = [1];
r = 2.4; %SDA's author suggested the value between 2.1~2.6

I = eye(2);
G = B*inv(R)*transpose(B);
H = transpose(C)*Q*C;
A_trans = transpose(A);
A_r = A - (r*I);

iteration_times = 0;

%solve CARE with SDA
A_hat_last = I + 2*r*inv(A_r + G*inv(transpose(A_r))*H);
G_hat_last = 2*r*inv(A_r)*G*inv(transpose(A_r) + H*inv(A_r)*G);
H_hat_last = 2*r*inv(transpose(A_r) + H*inv(A_r)*G)*H*inv(A_r);

while 1
    iteration_times = iteration_times + 1;
    
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

X_SDA = H_Hat_new;

%print SDA iteration times
disp(['iteration times = ' iteration_times]);
disp(iteration_times);

%print SDA result
disp("X (SDA):");
disp(X_SDA);

%solve CARE with MATLAB built-in function
[X_MATLAB, L_dummy, G_dummy] = care(A, B, H, R);

%print MATLAB CARE result
disp("X (MATLAB):");
disp(X_MATLAB);

%check if X really fits the CARE
CARE_SDA = transpose(A)*X_SDA + X_SDA*A - X_SDA*B*inv(R)*transpose(B)*X_SDA + transpose(C)*Q*C;
disp("Use X_SDA to calculate CARE:");
disp(CARE_SDA);
disp(norm(CARE_SDA));

%check if X really fits the CARE
CARE_MATLAB = transpose(A)*X_MATLAB + X_MATLAB*A - X_MATLAB*B*inv(R)*transpose(B)*X_MATLAB + transpose(C)*Q*C;
disp("Use X_MATLAB to calculate CARE:");
disp(CARE_MATLAB);
disp(norm(CARE_MATLAB));
