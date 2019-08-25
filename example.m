% Author: Haedong Kim
% Date: 2019-08-25
% It is an example of how to use HH_model function

%% run the model
% clear the environment
clear
close
clc

% set hyper-parameter I
I(1:500) = 50; 
I(501:2000) = 0; 
I(2001:10000) = 50;
step_size = 0.01;

% run the model
[V, g_Na, g_K] = HH_model(I, 0, step_size);

%% plot the results
x_lb = sprintf('Time Unit %s sec', step_size);

% V voltage
figure(1)
plot(V)
legend('Voltage')
ylabel('mV')
xlabel(x_lb)
title('Action Potential over Time')

% sodium conductance g_Na
figure(2)
plot(g_Na)
legend('Sodium Conductance')
xlabel(x_lb)
title('Sodium Conductance over Time')

% potassium conductance g_K
figure(3)
plot(g_K)
legend('Potassium Conductance')
xlabel(x_lb)
title('Potassium Conductance over Time')

% sodium & potassium conductances on the same plot
figure(4)
p1 = plot(g_Na);
hold on
p2 = plot(g_K);
legend([p1, p2], 'Sodium Conductance', 'Potassium Conductance')
xlabel(x_lb)
title('Ion Conductances over Time')
hold off
