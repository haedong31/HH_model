%% constants
clc
close all
clear variables

% conductances, cG_Kacitance
G_Na_max=120;
G_K_max=36;
G_L=0.3;
capa=1;

% resting potentials
V_Na=115;
V_K=-12;
V_L=10.613;

% initial values
V0=0;
n0=0;
m0=0;
h0=1;

%% run simulation as applying different external stimulation currents
time=0:0.01:1;
amplitudes=10:10:50;
sim_outs=cell(1,length(amplitudes));
I_m_container=zeros(length(time),length(amplitudes)+1);
I_m_container(:,1)=time;
for n=1:length(amplitudes)
    I_m=zeros(size(time));
    I_m=[time.' I_m.'];
    I_m(30:70, 2)=amplitudes(n);
    I_m_container(:,n+1)=I_m(:,2);
    
    sim_outs{n}=sim('HH_ap.slx');
end

%% plots
% plot action potential
figure(1)
for n=1:length(sim_outs)
    outputs=sim_outs{n}.get('yout');
    AP=outputs.getElement('AP');

    hold on
    plot(AP.Values.Time, AP.Values.Data, 'LineWidth', 2)
    hold off
end

legend('I_m = 10 uA/cm^2','I_m = 20 uA/cm^2','I_m = 30 uA/cm^2','I_m = 40 uA/cm^2','I_m = 50 uA/cm^2')
title('Action Potentioal')
xlabel('mSec')
ylabel('mV')

% plot I_Na
figure(2)
for n=1:length(sim_outs)
    outputs=sim_outs{n}.get('yout');
    I_Na=outputs.getElement('I_Na');
    
    hold on
    plot(I_Na.Values.Time, I_Na.Values.Data, 'LineWidth', 2)
    hold off
end

legend('I_m = 10 uA/cm^2','I_m = 20 uA/cm^2','I_m = 30 uA/cm^2','I_m = 40 uA/cm^2','I_m = 50 uA/cm^2')
title('Sodium Channel Current')
xlabel('mSec')
ylabel('uA/cm^2')

% plot I_K
figure(3)
for n=1:length(sim_outs)
    outputs=sim_outs{n}.get('yout');
    I_K=outputs.getElement('I_K');
    
    hold on
    plot(I_K.Values.Time, I_K.Values.Data, 'LineWidth', 2)
    hold off
end

legend('I_m = 10 uA/cm^2','I_m = 20 uA/cm^2','I_m = 30 uA/cm^2','I_m = 40 uA/cm^2','I_m = 50 uA/cm^2')
title('Potassium Channel Current')
xlabel('mSec')
ylabel('uA/cm^2')

% plot G_Na
figure(4)
for n=1:length(sim_outs)
    outputs=sim_outs{n}.get('yout');
    G_Na=outputs.getElement('G_Na');
    
    hold on
    plot(G_Na.Values.Time, G_Na.Values.Data, 'LineWidth', 2)
    hold off
end

legend('I_m = 10 uA/cm^2','I_m = 20 uA/cm^2','I_m = 30 uA/cm^2','I_m = 40 uA/cm^2','I_m = 50 uA/cm^2')
title('Sodium Channel Conductance')
xlabel('mSec')
ylabel('ms/cm')


%% plot G_K
figure(5)
for n=1:length(sim_outs)
    outputs=sim_outs{n}.get('yout');
    G_K=outputs.getElement('G_K');
    
    hold on
    plot(G_K.Values.Time, G_K.Values.Data, 'LineWidth', 2)
    hold off
end

legend('I_m = 10 uA/cm^2','I_m = 20 uA/cm^2','I_m = 30 uA/cm^2','I_m = 40 uA/cm^2','I_m = 50 uA/cm^2')
title('Potassium Channel Conductance')
xlabel('mSec')
ylabel('ms/cm')


%% plot I_m
figure(6)
for n=1:length(amplitudes)
    hold on
    plot(I_m_container(:,1),I_m_container(:,n+1),'LineWidth', 2)
    hold off
end

legend('I_m = 10 uA/cm^2','I_m = 20 uA/cm^2','I_m = 30 uA/cm^2','I_m = 40 uA/cm^2','I_m = 50 uA/cm^2')
title('Applied External Stimulation Current')
xlabel('mSec')
ylabel('uA/cm^2')
