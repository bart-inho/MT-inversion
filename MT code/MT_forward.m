% Magnetotelluric (MT) forward model

clear, close

% Parameters
load freq.mat;
freq = freq'; % [Hz] frequency of signal
T = 1./freq; % [s] period of signal
omega = 2*pi*freq; % [Hz] angular frequency of signal
z = [0 200 400 1000]'; % [m] depths
sigma = [1e-3 1e-1 1e-2 1e-5]'; % [S/m] electrical conductivites
fs = 13; % Fontsize

% Variables
mu_0 = 4*pi*1e-7; % [H/m] magnetic permeability of free space

[C] = Wait_recursion(omega,z,sigma,mu_0); % C-response [1/m]
rho_a = (real(C).^2+imag(C).^2)*mu_0.*omega; % Apparent resistivity [Ohm.m]
phi = atand(imag(C)/real(C)) + 90; % Impedance phase [deg]

% Regrouping data
% data = [freq';rho_a';phi'];

figure(1),clf
loglog(T,rho_a,'-o')
% plot(T, rho_a, '-o')
grid on
xlabel('Period T [s]','FontSize',fs)
ylabel('Apparent resistivity \rho_a [\Omega\cdotm]','FontSize',fs)
title('MT forward model with synthetic datas','FontSize',fs)




