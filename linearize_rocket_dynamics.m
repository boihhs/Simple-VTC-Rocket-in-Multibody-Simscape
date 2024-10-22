function [A, B, xChange] = linearize_rocket_dynamics(x0, u0)
    % Symbolic variables for state and input
    syms r_x r_y r_z v_x v_y v_z q0 q1 q2 q3 omega1 omega2 omega3 real
    syms fx fy fz real
    
    % State and input vectors
    x = [r_x; r_y; r_z; q0; q1; q2; q3; v_x; v_y; v_z; omega1; omega2; omega3];
    u = [fx; fy; fz];
    
    m = 4; % mass in Kg
    h = 1; % height in m
    g = 9.81; % gravity constant
    Radius = 95; % Radius in m
    Radius = Radius * 0.001;
    %Interia
    J = [ [m/4*(1/3*h^2+Radius^2) 0 0]; [0 m/4*(1/3*h^2+Radius^2) 0]; [0 0 (m/2*Radius^2)]];

    H = [0, 0, 0;
        eye(3)];

    r_offset = [0; 0; -h/2];


    r = [r_x; r_y; r_z];
    q = [q0; q1; q2; q3];
    v = [v_x; v_y; v_z];
    w = [omega1; omega2; omega3];

    goalq = [1; 0; 0; 0];
    changeGoalq = inverseCayleyMap(L(goalq)'*q);
    xChange_k = [r; changeGoalq; v; w];



    dr = v;
    dq = .5*L(q)*H*w;
    dv = 1/m * (rotateVector(q, u)) + [0; 0; -9.81];
    dw = J^-1 * (cross(r_offset, u) - cross(w, J*w));
    
   
    
    % State equations
    f = [dr; dq; dv; dw];
    
    % Jacobians (A = df/dx, B = df/du)
    A_sym = jacobian(f, x);
    B_sym = jacobian(f, u);
   
    
    % Evaluate E(x) Jacobian transformation matrix
    E_x = E_matrix(q);          % E(x) function to handle quaternion Jacobians
    
    A_k = E_x' * A_sym * E_x;
    B_k = E_x' * B_sym;
    
    % Substitute the operating point values
    A = double(subs(A_k, [x; u], [x0; u0]));
    B = double(subs(B_k, [x; u], [x0; u0]));
    
    xChange = double(subs(xChange_k, [x; u], [x0; u0]));
    
  
end

function L = L(q)
    qs = q(1);
    qv = q(2:4);
    L = [qs, -qv';
        qv, qs*eye(3) + skew(qv)];
end

function R = R(q)
    qs = q(1);
    qv = q(2:4);
    R = [qs, -qv';
        qv, qs*eye(3) - skew(qv)];
end

function A = rotateVector(q, r)
    H = [0, 0, 0;
        eye(3)];
    A = H'*L(q)*R(q)'*H*r;
end

function x = skew(q)
    x = [0, -q(3), q(2);
        q(3), 0, -q(1);
        -q(2), q(1), 0];
end

function G = G(q)
    qs = q(1);
    qv = q(2:4);
    G = [-qv';
        qs*eye(3) + skew(qv)];
end

function E = E_matrix(q)
    
    E = [eye(3), zeros(3, 9);
        zeros(4, 3), G(q), zeros(4, 6);
        zeros(3, 6), eye(3), zeros(3, 3);
        zeros(3, 9), eye(3)];
end

function phi = inverseCayleyMap(q)
    qs = q(1);
    qv = q(2:4);
    phi = qv/qs;
end

