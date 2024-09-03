% Run generateData.m to generate required files for this code

load('potential_Rtilde_c_1.mat')
potential_c1 = potential_Rtilde;
clear potential_Rtilde

load('potential_Rtilde_c_2.mat')
potential_c2 = potential_Rtilde;
clear potential_Rtilde

load('potential_Rtilde_c_3.mat')
potential_c3 = potential_Rtilde;
clear potential_Rtilde

load('potential_Rtilde_c_4.mat')
potential_c4 = potential_Rtilde;
clear potential_Rtilde

load('potential_Rtilde_c_5.mat')
potential_c5 = potential_Rtilde;
clear potential_Rtilde

load('potential_Rtilde_c_6.mat')
potential_c6 = potential_Rtilde;
clear potential_Rtilde

load('time.mat')

%%

lw = 1.1;

figure()
plot(time, potential_c1, 'Linewidth', lw)
hold on
plot(time, potential_c2, 'Linewidth', lw)
hold on
plot(time, potential_c3, 'Linewidth', lw)
hold on
plot(time, potential_c4, 'Linewidth', lw)
hold on
plot(time, potential_c5, 'Linewidth', lw)
hold on
plot(time, potential_c6, 'Linewidth', lw)

legend('$c_1 = 0.84$', '$c_1 = 0.86$', '$c_1 = 0.88$', '$c_1 = 0.90$', '$c_1 = 0.92$', '$c_1 = 0.94$', 'Interpreter', 'Latex')
xlabel("$t \: [s]$", 'Interpreter', 'latex')
ylabel({'$|\tilde{R}|^2_I $'}, 'interpreter', 'latex')
ax = gca;
ax.FontSize = 20;
grid on
pbaspect([3 1 1])

axes('position',[.4 .4 .25 .25])
box on % put box around new pair of axes
indexOfInterest = (time < 2.9) & (time > 2.8); % range of t near perturbation
plot(time(indexOfInterest),potential_c1(indexOfInterest), 'Linewidth', lw) % plot on new axes
hold on
plot(time(indexOfInterest),potential_c2(indexOfInterest), 'Linewidth', lw) % plot on new axes
hold on
plot(time(indexOfInterest),potential_c3(indexOfInterest), 'Linewidth', lw) % plot on new axes
hold on
plot(time(indexOfInterest),potential_c4(indexOfInterest), 'Linewidth', lw) % plot on new axes
hold on
plot(time(indexOfInterest),potential_c5(indexOfInterest), 'Linewidth', lw) % plot on new axes
hold on
plot(time(indexOfInterest),potential_c6(indexOfInterest), 'Linewidth', lw) % plot on new axes
grid on
axis tight
