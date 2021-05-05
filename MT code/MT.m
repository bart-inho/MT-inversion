% MT

clear

load freq.mat

% Parameters
f = freq;  % [Hz] frequency of signal
omega = 2*pi.*f; % [Hz] angular frequency of signal
z = linspace(0,200,length(f))'; % [m] depth of interfaces between layers
sigma = 2e-2*ones(length(f),1); % [S/m] electrical conductivity of layers

% Variables
mu_0 = 4*pi*1e-7; % [H/m] magnetic permeability

% Wait's recursion formula
[C,rho_a] = Wait_recursion(omega,z,sigma,mu_0);

figure(1),clf
plot(rho_a,omega,'-o')
xlabel('\rho_a')
ylabel('\omega')



