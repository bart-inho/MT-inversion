% Magnetotelluric (MT) forward model
% N-layered half-space

clear

% Parameters
load freq.mat;
freq = freq'; % [Hz] frequency of signal
T = 1./freq; % [s] period of signal
omega = 2*pi*freq; % [Hz] angular frequency of signal
z = [0 5000 10000]'; % [m] depths
sigma = [0.01 3 0.01]'; % [S/m] electrical conductivites
fs = 13; % Fontsize

% Variables
mu_0 = 4*pi*1e-7; % [H/m] magnetic permeability of free space

[C] = Wait_recursion(omega,z,sigma,mu_0); % C-response [m]
rho_a = abs(C).^2*mu_0.*omega; % Apparent resistivity [Ohm.m]
phi = atand(imag(C)/real(C)) + 90; % Impedance phase [deg]

% Regrouping data
% data = [freq';rho_a';phi'];

figure
loglog(T,rho_a,'.-')
% plot(T, rho_a, '-o')
grid on
xlabel('Period T [s]','FontSize',fs)
ylabel('Apparent resistivity \rho_a [\Omega\cdotm]','FontSize',fs)
title('MT forward model with synthetic data','FontSize',fs)



