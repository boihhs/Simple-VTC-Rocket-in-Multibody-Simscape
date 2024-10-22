function U_opt = lqrC(xin, uin, t)

Ts = .05;
m = 4; % mass in Kg
horizon = 100;
h = 1; % height in m
g = 9.81; % gravity constant
Radius = 95; % Radius in m
Radius = Radius * 0.001;
moi_tensor = [ [m/4*(1/3*h^2+Radius^2) 0 0]; [0 m/4*(1/3*h^2+Radius^2) 0]; [0 0 (m/2*Radius^2)]];
rbody = [0, 0, -h/2]';
nx = 12;
nu = 3;
[Ac, Bc, xChange] = linearize_rocket_dynamics(xin, uin);

Q = eye(nx) * 1000;

R = eye(nu);

% Correct Augmented Matrix for Discretization (16x16)
augmented_matrix = [Ac, Bc; zeros(3, 12 + 3)] * Ts;

% Compute matrix exponential
exp_augmented = expm(augmented_matrix);

% Extract Ad and Bd
Ad = exp_augmented(1:12, 1:12);
Bd = exp_augmented(1:12, 13:15);


K = zeros(size(Bd, 2), size(Ad, 1), horizon);    
S = Q;
disp(Ad*Bd);
disp(Ad);
for i = horizon:-1:1
    %disp((Bd' * S*Ad))
    K(:, :, i) = (Bd' * S * Bd + R)^-1*(Bd' * S * Ad);
    
    S = Q + Ad'*S*Ad - (Ad'*S*Bd)*(R+ Bd'*S*Bd)^-1*(Bd'*S*Ad);
end
%disp(K(:, :, 2));
%disp(xChange);
U_opt = -K(:, :, 1)*xChange;

disp("Time: " + t)
disp(U_opt);



end
