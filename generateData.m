% run this file before executing compareConvergenceRates.m

start_time = 0.;
end_time = 10.;
N = 1000*(end_time - start_time);
time = linspace(start_time, end_time, N);
h = time(2)-time(1);
save('time.mat','time')

%% Initializations
axis = [1;0;0]; % fix this axis of rotation for all elements of SO(3), resulting in SO(2).
init_Rtrue = fun_axisangle((pi/180.) * 0.  , axis); % initial value of true rotation matrix
init_Rhat  = fun_axisangle((pi/180.) * 180,  axis); % initial value of attitude estimate
init_q = 0;


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

q = zeros(N,1);
q(1) = init_q;

Omega_true = zeros(3,N); % true angular velocity 
for i=1:1:N
    Omega_true(:,i) = sin(time(i))*axis;
end

Omega_y = Omega_true; % measurement angular velocity

jumps = zeros(N,1); % number of jumps

%% Filter parameters

kp = 1.0; % Observer gain
c0 = fun_potential(fun_axisangle((pi/180.)*150.0, axis)); % actually c0^2
% c1 = fun_potential(fun_axisangle((pi/180.)*120.0, axis)); % actually c1^2

c1_array = [0.84^2; 0.86^2; 0.88^2; 0.9^2; 0.92^2; 0.94^2];

% Rstar_angle = (pi/180.)*90.; 
% norm_Rstar = 0.5*(1 - cos(Rstar_angle)); % actually |R^*|^2
norm_Rstar = (0.58)^2; % actually |R^*|^2
Rstar_angle = acos(1 - 2*norm_Rstar);
Rstar = fun_axisangle(Rstar_angle, axis);


for r = 1:1:6
if (norm_Rstar < c1_array(r)) || (norm_Rstar < c0) || (c1_array(r) < c0) || (c1_array(r) + norm_Rstar > 1)
    if (4/c1_array(r)) < 3 + 1/norm_Rstar
        continue
    else
        warning('No exponential bound');
    end
else
    error('Parameters not set appropriately.');
end
end


%% Main loop
% Store |\tilde R|^2 for the hybrid filter

potential_Rtilde = zeros(N,1);
potential_Rtilde(1) = fun_potential(Rtilde(:,:,1));

% Store |\tilde R|^2 for PCF
potential_Rtilde_nonHybrid = zeros(N,1);
potential_Rtilde_nonHybrid(1) = fun_potential(Rtilde_nonHybrid(:,:,1));

Vdot = zeros(N,1);

jump_times = [];
jump_index = [];
noise = 1; % 1 == noisy measurements, 0 == no noise

for p=1:1:length(c1_array) % each loop for a value of c1
    c1 = c1_array(p);
    for i=1:1:(N-1)
        Rtrue(:,:,i+1) = fun_rotationPropagation(Rtrue(:,:,i), Omega_true(:,i), h); % true kinematics propagation
    
        if noise == 1  
            Rtilde_nonHybrid_angle = real((acos(1 - 2*(potential_Rtilde_nonHybrid(i)))))*(180/pi);
            if potential_Rtilde_nonHybrid(i) < 0.999
                noise_omega = -1*sign(sin(Rtilde_nonHybrid_angle*pi/180.))*axis;
            else
                noise_omega = sin(Rtilde_nonHybrid_angle*pi/180.)*axis;
            end
            noise_angle = sin(2*pi*time(i)/1)* 10* (pi/180.);
            noise_rotm = fun_axisangle(noise_angle, axis);

        % if noise == 1
        %     noise_omega = 0.1*rand(3,1);
        %     noise_angle = -pi/100 + rand()*pi/50;
        %     noise_rotm = fun_axisangle(noise_angle, axis);
    
        elseif noise == 0
            noise_omega = zeros(3,1);
            noise_rotm = eye(3);
        else
            error("wrong input for noise")
        end

        Omega_y(:,i) = Omega_true(:,i) + noise_omega; % noisy omega
        Ry(:,:,i) = Rtrue(:,:,i)*noise_rotm; % noisy rotation
    
        % update estimates for the hybrid filter
        [Rhat(:,:,i+1), q(i+1), jumps(i+1,1)] = ...
            fun_hybridPCF(  q(i), Rhat(:,:,i), Ry(:,:,i), norm_Rstar,...
                            Omega_y(:,i), kp, c0, c1, h, jumps(i,1)   );

        % update estimates for PCF
        Rhat_nonHybrid(:,:,i+1) = fun_passiveComplementaryFilter(Rhat_nonHybrid(:,:,i),...
            Ry(:,:,i), Omega_y(:,i), kp, h);
                    
        % update filter error
        Rtilde(:,:,i+1) = Rhat(:,:,i+1)' * Rtrue(:,:,i+1);
        Rtilde_nonHybrid(:,:,i+1) = Rhat_nonHybrid(:,:,i+1)' * Rtrue(:,:,i+1);

        % calculate |\tilde R|^2 for the hybrid filter and PCF
        potential_Rtilde(i+1) = fun_potential(Rtilde(:,:,i+1));
        potential_Rtilde_nonHybrid(i+1) = fun_potential(Rtilde_nonHybrid(:,:,i+1));
    
        % calculate jump times
        if jumps(i+1) - jumps(i) ~= 0
            jump_times(end+1) = time(i);
            jump_index(end+1) = i;
        end
    end
    
    save(['potential_Rtilde_c_' num2str(p) '.mat'],'potential_Rtilde');
end