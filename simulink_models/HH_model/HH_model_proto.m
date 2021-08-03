%% simulator based on Hodgkin-Huxley model
function [V_m, g_Na, g_K] = HH_model_proto(I_m, V_0)
%{
This function retruns a simulation results of potential associated to
memebrane mechanism specified by Hodgkin-Huxley model.
:param I_m:
:return V_m:
:return g_Na:
:return_g_K:
%}

% constants
E_Na = -115; 
E_K = 12; 
E_L = -10.613;
G_Na = 120;  
G_K = 36; 
G_L = 0.3;
C_m = 1;

% initialization
V = V_0;

% solve a system of differential equations
syms m(t) n(t) h(t);
ode_m = diff(m, t) == alpha_m(V)*(1-m) - beta_m(V)*m;
ode_n = diff(n, t) == alpha_n(V)*(1-n) - beta_n(V)*n;
ode_h = diff(h, t) == alpha_h(V)*(1-h) - beta_h(V)*h;
odes = [ode_m; ode_n; ode_h];
dsolve(odes)
end


%% functions for gating variables
function am = alpha_m(V)
%{
:param V: diviation of the membrane voltage from the resting voltage
:return: transfer rate coefficient for m-particles from open to close
%}
am = (0.1*(25-V))/(exp((25-V)/10)-1);
end

function bm = beta_m(V)
%{
:param V: diviation of the membrane voltage from the resting voltage
:return: transfer rate coefficient for m-particles from close to open
%}
bm = 4/(exp(V/18));
end

function ah = alpha_h(V)
%{
:param V: diviation of the membrane voltage from the resting voltage
:return:
%}
ah = 0.07/exp(V/20);
end

function bh = beta_h(V)
%{
:param V: diviation of the membrane voltage from the resting voltage
:return:
%}
bh = 1/(exp((30-V)/10)+1);
end

function an = alpha_n(V)
%{
:param V: diviation of the membrane voltage from the resting voltage
:return:
%}
an = (0.01*(10-V))/(exp((10-V)/10)-1);
end

function bn = beta_n(V)
%{
:param V: diviation of the membrane voltage from the resting voltage
:return:
%}
bn = 0.125/exp(V/80);
end
