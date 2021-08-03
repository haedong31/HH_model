%% Hodgkin-Huxley simulator
function [V, g_Na, g_K] = HH_model(I, V_base, step_size)
%{
    This function retruns a simulation results of potential associated to
    memebrane mechanism specified by Hodgkin-Huxley model. Solutions are
    obtained by Euler's method.
    :param I: vector of total currents on axon
    :param V_base: base voltage
    :param step_size: step size of update in Euler's method
    :return V: vector of action potential
    :return g_Na: vector of sodium conductances
    :return_g_K: vecotr of potassium conductances
%}

% number of elements of vector I
num_steps = numel(I);

% constants
E_Na = 115;
E_K = -12;
E_L = 10.613;
G_Na = 120;
G_K = 36;
G_L = 0.3;
C_m = 1;

% set empty vectors that will contain values of interest
V = zeros(1, num_steps);
m = zeros(1, num_steps);
n = zeros(1, num_steps);
h = zeros(1, num_steps);
alpha_m = zeros(1, num_steps);
beta_m = zeros(1, num_steps);
alpha_n = zeros(1, num_steps);
beta_n = zeros(1, num_steps);
alpha_h = zeros(1, num_steps);
beta_h = zeros(1, num_steps);
I_Na = zeros(1, num_steps);
I_K = zeros(1, num_steps);
I_L = zeros(1, num_steps);
I_m = zeros(1, num_steps);
g_Na = zeros(1, num_steps);
g_K = zeros(1, num_steps);

% transfer coefficient rates in resting state
alpha_n0 = .01 * ( (10-V_base) / (exp((10-V_base)/10)-1) );
beta_n0 = .125*exp(-V_base/80);
alpha_m0 = .1*( (25-V_base) / (exp((25-V_base)/10)-1) );
beta_m0 = 4*exp(-V_base/18);
alpha_h0 = .07*exp(-V_base/20);
beta_h0 = 1/(exp((30-V_base)/10)+1);

% initial values of V, m, n, and h in resting state
V(1) = V_base;
m(1) = alpha_m0/(alpha_m0+beta_m0);
n(1) = alpha_n0/(alpha_n0+beta_n0);
h(1) = alpha_h0/(alpha_h0+beta_h0);

for i=1:num_steps-1
    % evaluate transfer coefficients
    alpha_n(i) = 0.01*((10-V(i)) / (exp((10-V(i))/10)-1));
    beta_n(i) = 0.125*exp(-V(i)/80);
    alpha_m(i) = 0.1*((25-V(i)) / (exp((25-V(i))/10)-1));
    beta_m(i) = 4*exp(-V(i)/18);
    alpha_h(i) = 0.07*exp(-V(i)/20);
    beta_h(i) = 1 / (exp((30-V(i))/10)+1);
    
    % calculate currents
    I_Na(i) = G_Na * m(i)^3 * h(i) * (V(i)-E_Na);
    I_K(i) = G_K * n(i)^4 * (V(i)-E_K);
    I_L(i) = G_L * (V(i)-E_L);
    I_m(i) = I(i) - (I_Na(i) + I_K(i) + I_L(i));
    
    % estimate the next values
    V(i+1) = V(i) + step_size/C_m*(I_m(i));
    m(i+1) = m(i) + step_size*(alpha_m(i)*(1-m(i)) - beta_m(i)*m(i));
    n(i+1) = n(i) + step_size*(alpha_n(i)*(1-n(i)) - beta_n(i)*n(i));
    h(i+1) = h(i) + step_size*(alpha_h(i)*(1-h(i)) - beta_h(i)*h(i));
end

% transfer rate coefficients at the last time period
alpha_m(num_steps) = (0.1*(25-V(num_steps)))/(exp((25-V(num_steps))/10)-1);
beta_m(num_steps) = 4/(exp(V(num_steps)/18));
alpha_h(num_steps) = 0.07/exp(V(num_steps)/20);
beta_h(num_steps) = 1/(exp((30-V(num_steps))/10)+1);
alpha_n(num_steps) = (0.01*(10-V(num_steps)))/(exp((10-V(num_steps))/10)-1);
beta_n(num_steps) = 0.125/exp(V(num_steps)/80);

% currents at the last time period
I_Na(num_steps) = G_Na*m(num_steps)^3*h(num_steps)*(V(num_steps)-E_Na);
I_K(num_steps) = G_K*n(num_steps)^4*(V(num_steps)-E_K);
I_L(num_steps) = G_L*(V(num_steps)-E_L);
I_m(num_steps) = I(num_steps) - (I_Na(num_steps) + I_K(num_steps) + I_L(num_steps));

% calculate conductances
for i = 1:num_steps
    g_Na(i) = I_Na(i)/(V(i) - E_Na);
    g_K(i) = I_K(i)/(V(i) - E_K);
end

end
