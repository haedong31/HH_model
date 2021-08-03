%% example 1 Single ODE
syms y(t);
ode = diff(y, t) == t*y;
ySol(t) = dsolve(ode);
disp(ySol(t))

%% example 2 System of differential equations
syms u(t) v(t);
ode1 = diff(u, t) == 3*u + 4*v;
ode2 = diff(v, t) == -4*u + 3*v;
odes = [ode1; ode2];
sols = dsolve(odes);
disp(sols)

%% hh
% constants (parameters)
E_Na = -115; 
E_K = 12; 
E_L = -10.613;
G_Na = 120;  
G_K = 36; 
G_L = 0.3;
C_m = 1;

% initial or given values
I = 50;

alpha_m = (0.1*(25-V))/(exp((25-V)/10)-1);
beta_m = 4/(exp(V/18));
alpha_h = 0.07/exp(V/20);
beta_h = 1/(exp((30-V)/10)+1);
alpha_n = (0.01*(10-V))/(exp((10-V)/10)-1);
beta_n = 0.125/exp(V/80);

syms V(t) m(t) n(t) h(t);
ode_V = diff(V, t) == (I - (G_K*n^4*(V-E_K) + G_Na*m^3*h*(V-E_Na) + G_L*(V-E_L))) / C_m;
ode_m = diff(m, t) == alpha_m(V)*(1-m) - beta_m(V)*m;
ode_n = diff(n, t) == alpha_n(V)*(1-n) - beta_n(V)*n;
ode_h = diff(h, t) == alpha_h(V)*(1-h) - beta_h(V)*h;
odes = [ode_V; ode_m; ode_n; ode_h];

% cannot find explicit solution
dsolve(odes)
