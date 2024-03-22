% Hybrid passive complementary filter implementation in SO(3)
% Comment out appropriate lines of code to run Case 1 and Case 2. 

clear;
clc;

start_time = 0.;
end_time = 10.;
N = 1000*(end_time - start_time);
time = linspace(start_time, end_time, N);
h = time(2)-time(1);

%% Initializations
axis = [1;0;0];
init_Rtrue = fun_axisangle((pi/180.) * 0.  ,  [1;0;0]); % initial value of true rotation matrix
% init_Rhat  = fun_axisangle((pi/180.) * 120.,  [0;1;0]); % Case 1: attitude estimate initial condition
init_Rhat  = fun_axisangle((pi/180.) * 180.,  [0;1;0]); % Case 2: attitude estimate initial condition for 


Rtrue = zeros(3,3,N); % true rotation matrix
Rtrue(:,:,1) = init_Rtrue;

Ry = zeros(3,3,N); % measured rotation matrix
Ry(:,:,1) = Rtrue(:,:,1);

Rhat = zeros(3,3,N); % estimated rotation matrix - hybrid filter
Rhat(:,:,1) = init_Rhat;

Rhat_nonHybrid = zeros(3,3,N); % estimated rotation matrix - PCF
Rhat_nonHybrid(:,:,1) = init_Rhat;

Rtilde = zeros(3,3,N); % filter error - hybrid filter
Rtilde(:,:,1) = Rhat(:,:,1)' * Rtrue(:,:,1);

Rtilde_nonHybrid = zeros(3,3,N); % filter error - PCF
Rtilde_nonHybrid(:,:,1) = Rhat_nonHybrid(:,:,1)' * Rtrue(:,:,1);

Rtilde_y = zeros(3,3,N); % attitude measurement error
Rtilde_y(:,:,1) = Rtilde(:,:,1);

Omega_true = zeros(3,N); % true angular velocity 
for i=1:1:N
    Omega_true(1,i) = sin(time(i));
    Omega_true(2,i) = cos(time(i));
    Omega_true(3,i) = 1.;
end

Omega_y = Omega_true; % measurement angular velocity

jumps = zeros(N,1); % number of jumps

%% Filter parameters

kp = 1.0; % Observer gain
c = fun_potential(fun_axisangle((pi/180.)*150.0, axis)); % actually c^2

%% Main loop

potential_Rtilde = zeros(N,1);
potential_Rtilde(1) = fun_potential(Rtilde(:,:,1));

potential_Rtilde_nonHybrid = zeros(N,1);
potential_Rtilde_nonHybrid(1) = fun_potential(Rtilde_nonHybrid(:,:,1));

jump_times = [];
jump_index = [];
noise = 1; % 1 == noisy, 0 == no noise

for i=1:1:N-1
    Rtrue(:,:,i+1) = fun_rotationPropagation(Rtrue(:,:,i), Omega_true(:,i), h); % true kinematics propagation
    
    if noise == 1  

%         % Case 1
%         noise_omega_axis = rand(3,1);
%         noise_omega_axis = noise_omega_axis/norm(noise_omega_axis); 
%         max_omega = 0.5;
%         noise_omega = max_omega * rand() * noise_omega_axis;
% 
%         max_noise_angle = (pi/180.)*10.;
%         noise_angle = -max_noise_angle + 2*max_noise_angle*rand();
%         noise_axis = rand(3,1);
%         noise_axis = noise_axis/norm(noise_axis);
%         noise_rotm = fun_axisangle(noise_angle, axis);

%         % Case 2
        omega_Rtilde = fun_vex(fun_skewSymmetricMatrix(Rtilde_nonHybrid(:,:,i)));
        noise_omega = -sign(omega_Rtilde);
        noise_angle = sin(2*pi*time(i)/1)* 10* (pi/180.);
        noise_rotm = fun_axisangle(noise_angle, axis);

    elseif noise == 0
        noise_omega = zeros(3,1);
        noise_rotm = eye(3);
    else
        error("wrong input for noise")
    end
    
    Omega_y(:,i) = Omega_true(:,i) + noise_omega; % noisy omega
    Ry(:,:,i) = Rtrue(:,:,i)*noise_rotm; % noisy rotation

    % update estimates
    [Rhat(:,:,i+1), jumps(i+1,1)] = ...
    fun_hybridPCF( Rhat(:,:,i), Ry(:,:,i), Omega_y(:,i), kp, c, h, jumps(i,1) );

    Rhat_nonHybrid(:,:,i+1) = fun_passiveComplementaryFilter(Rhat_nonHybrid(:,:,i),...
        Ry(:,:,i), Omega_y(:,i), kp, h);

    Rtilde(:,:,i+1) = Rhat(:,:,i+1)' * Rtrue(:,:,i+1);
    Rtilde_nonHybrid(:,:,i+1) = Rhat_nonHybrid(:,:,i+1)' * Rtrue(:,:,i+1);

    % calculate potential_Rtilde
    potential_Rtilde(i+1) = fun_potential(Rtilde(:,:,i+1));
    potential_Rtilde_nonHybrid(i+1) = fun_potential(Rtilde_nonHybrid(:,:,i+1));
    
    % calculate jump times
    if jumps(i+1) - jumps(i) ~= 0
        jump_times(end+1) = time(i);
        jump_index(end+1) = i;
    end
end

number_of_jumps = length(jump_index);


%% Plots for Case 1
clf;

% figure(1)
% plot(time(2:end), potential_Rtilde(2:end,:), 'LineWidth', 2)
% xlabel("$t \: [s]$", 'Interpreter', 'latex')
% ylabel({'$|\tilde{R}|^2_I $'}, 'interpreter', 'latex')
% ax = gca;
% ax.FontSize = 20;
% grid on
% legend('Hybrid Filter on $\textrm{SO}(3)$', 'Interpreter', 'latex')
% xlabel("$t \: [s]$", 'Interpreter', 'latex')
% ylabel({'$|\tilde{R}|^2_I $'}, 'interpreter', 'latex')
% ax = gca;
% ax.FontSize = 20;
% grid on
% if number_of_jumps == 0
%     return;
% else
%     for i=1:1:number_of_jumps
%         ind = jump_index(i);
%         hold on;
%         plot(time(ind), potential_Rtilde(ind), 'color', 'red', 'marker','x', 'linewidth', 2, 'MarkerSize',12);
%         hold on;
%         plot(time(ind+1), potential_Rtilde(ind+1), 'color', 'red', 'marker','o', 'linewidth', 1, 'MarkerSize',9);
%         hold on;
%         plot([time(ind) time(ind+1)], [potential_Rtilde(ind) potential_Rtilde(ind+1)], '--r', 'LineWidth',2);
%     end
% end


%% Plots for Case 2

figure(1)
plot(time, potential_Rtilde_nonHybrid, 'color', 'black', 'LineWidth', 2);
hold on;
plot(time(2:end), potential_Rtilde(2:end,:),'color', 'blue', 'LineWidth', 2)
xlabel("$t \: [s]$", 'Interpreter', 'latex')
ylabel({'$|\tilde{R}|^2_I $'}, 'interpreter', 'latex')
ax = gca;
ax.FontSize = 20;
grid on
if number_of_jumps == 0
    return;
else
    for i=1:1:number_of_jumps
        ind = jump_index(i);
        hold on;
        plot(time(ind), potential_Rtilde(ind), 'color', 'red', 'marker','x', 'linewidth', 2, 'MarkerSize',12);
        hold on;
        plot(time(ind+1), potential_Rtilde(ind+1), 'color', 'red', 'marker','o', 'linewidth', 1, 'MarkerSize',9);
        hold on;
        plot([time(ind) time(ind+1)], [potential_Rtilde(ind) potential_Rtilde(ind+1)], '--r', 'LineWidth',2);
    end
end
hold on;
plot(time, potential_Rtilde_nonHybrid, 'color', 'black', 'LineWidth', 2);
ylim([0 1.1])
legend('Passive Complementary Filter', 'Hybrid Filter on $\textrm{SO}(3)$', '', '', '','', 'Interpreter', 'latex')

xlabel("$t \: [s]$", 'Interpreter', 'latex')
ylabel({'$|\tilde{R}|^2_I $'}, 'interpreter', 'latex')
ax = gca;
ax.FontSize = 20;
grid on
