h = 0.0001;

N = 10000;
axis = [1;0;0];

time_f = 0;
time_g = 0;

for i=1:1:N
    theta_hat = 2*pi*rand();
    theta_meas = 2*pi*rand();

    kp = rand();
    norm_Rstar = 0.1 + 0.4*rand();
    
    Rhat = fun_axisangle(theta_hat, axis);
    Ry = fun_axisangle(theta_meas, axis);
    Omega_y = rand(3,1);

    c1 = rand();
    c0 = c1 + (1-c1)*rand();

    q = randi([0,1]);

    f = @() fun_passiveComplementaryFilter(Rhat, Ry, Omega_y, kp, h);
    g = @() fun_hybridPCF(q, Rhat, Ry, norm_Rstar, Omega_y, kp, c0, c1, h);

    time_f = time_f + timeit(f);
    time_g = time_g + timeit(g);
end

time_f = time_f/N
time_g = time_g/N
