%% 1
% TODO: add frames sym(pi)cz

theta = sym('theta', [3 1]);
syms a b c d e f g;

% ans
% columns: theta d a alpha
dh = [0, d + b, -c, -sym(pi)/2; % FR->F0
    theta(1) - sym(pi)/2, 0, e, 0; % F1
    theta(2), 0, f, 0; % F2
    theta(3), 0, g, 0; % F3
    -sym(pi)/2, 0, 0, -sym(pi)/2] % FT

%% 2

T = sym(zeros(4,4,size(dh, 1)));
for joint=1:size(dh, 1)
    T(:,:,joint) = dh2mat(dh(joint,1),dh(joint,2),dh(joint,3),dh(joint,4));
end

% ans
T_R_T = T(:,:,1)*T(:,:,2)*T(:,:,3)*T(:,:,4)*T(:,:,5)

%% 3
x_t = T_R_T * [0;0;0;1];
x_t = x_t(1:3,:);
J_upper = jacobian(x_t, theta);
J_lower = [0,0,0;
           1,1,1; % thetas contribute to rotation in y
           0,0,0];
% ans
J = vertcat(J_upper, J_lower)
% ans
J_reduced_dof = J([1 3 5], :)

%% 4a
T_R_T_numeric = subs(T_R_T, [b, c, d, e, f, g], ...
                [361, 250, 380, 328, 323, 82.4]);
T_R_T_numeric = subs(T_R_T_numeric, theta, [sym(pi)/3; sym(pi)/2; sym(pi)/3]);

% ans
double(T_R_T_numeric)

%% 4b

J_numeric = subs(J, [b, c, d, e, f, g], ...
                [361, 250, 380, 328, 323, 82.4]); % mm
J_numeric = subs(J_numeric, theta, [sym(pi)/3; sym(pi)/2; sym(pi)/3]); % rad

q_dot_numeric = [sym(pi)/4; sym(pi)/4; sym(pi)/4]; % rad/s
x_dot_numeric = J_numeric * q_dot_numeric; % mm, rad/s

x_dot_numeric_deg = double(x_dot_numeric);
x_dot_numeric_deg(4:6,:) = rad2deg(x_dot_numeric_deg(4:6,:));

% ans
x_dot_numeric_deg % mm, deg/s

%% 4c
F_numeric = [50; 0; 0; 0; 0; 0]; % N, N*mm
joint_torque = J_numeric'*F_numeric;

% ans
double(joint_torque) % N*mm


%% 
syms theta_base r;
% xi_dot = sym('xi_dot', [3 1]); % x_dot y_dot theta_dot
phi_dot = sym('phi_dot', [4,1]); %left right omni_big omni_small


alpha = [sym(pi)/2; -sym(pi)/2; sym(pi)];
beta = [0; sym(pi); -sym(pi)/2];
l = [0.5*a; 0.5*a; c];

% equations stolen from lecture slides
% left hand side vectors
R_rob_0 = @(th) [cos(th), sin(th), 0; -sin(th), cos(th), 0; 0, 0, 1];
cons_fixed_roll = @(a, b, l) [sin(a+b), -cos(a+b), -l*cos(b)];
cons_fixed_slide = @(a, b, l) [cos(a+b), sin(a+b), -l*sin(b)];
cons_omni_roll = @(a, b, l) [sin(a+b), -cos(a+b), -l*cos(b)];
cons_omni_slide = @(a, b, l) [cos(a+b), sin(a+b), l*sin(b)];
J1_rolling = ...
    [cons_fixed_roll(alpha(1), beta(1), l(1)); % left roll
    cons_fixed_roll(alpha(2), beta(2), l(2)); % right roll
    cons_omni_roll(alpha(3), beta(3), l(3))] % omni roll
C1_sliding = ...
    [cons_fixed_slide(alpha(1), beta(1), l(1)); % left slide
    cons_fixed_slide(alpha(2), beta(2), l(2)); % right slide
    cons_omni_slide(alpha(3), beta(3), l(3))] % omni slide

J2_rolling = ...
    [r;
    r;
    r]
C2_sliding = ...
    [0;
    0;
    r*phi_dot(4)]

%% 6
rank([J1_rolling;C1_sliding])

%% 7
left_hand_side = [J1_rolling;C1_sliding]
right_hand_side = [J2_rolling.*phi_dot(1:3);C2_sliding]

% remove non-linearly-independent or non controllable equations
% We are not removing them. We just don't need them -Chan
left_hand_side = left_hand_side([1,2,4],:)
right_hand_side = right_hand_side([1,2,4],:)


% xi_0_dot = R_rob_0(theta_base)\left_hand_side\right_hand_side
xi_0_dot = inv(R_rob_0(theta_base))*inv(left_hand_side)*right_hand_side

%% 8
xi_0_dot_numeric = subs(xi_0_dot, [a, r, theta_base], [507, 143, sym(pi)/4]);
xi_0_dot_numeric = subs(xi_0_dot_numeric, phi_dot, [2*sym(pi);4*sym(pi);0;0]);
% ans
double(xi_0_dot_numeric) % mm/s rad/s














    