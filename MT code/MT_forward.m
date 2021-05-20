% Magnetotelluric (MT) forward model
% N-layered half-space
% 
% Camp de Geophysique d'Exploration
% Projet 5: Magnetotellurique
% Barthelemy Anhorn & Bruno Galissard de Marignac
% 
% First exercise

clear

% Parameters
load freq.mat;
freq = freq'; % [Hz] frequency of signal
T = 1./freq; % [s] period of signal
omega = 2*pi*freq; % [Hz] angular frequency of signal
z = [0 5000 10000]'; % [m] depths
thick = diff(z); % [m] thicknesses of layers, Note: dim(thick) = dim(z)-1
sigma = [0.01 3 0.01]'; % [S/m] electrical conductivities
rho = 1./sigma; % [Ohm.m] resistivities 
fs = 13; % Fontsize

% Variables
mu_0 = 4*pi*1e-7; % [H/m] magnetic permeability of free space


[C,rho_a,phi] = Wait_recursion(T,thick,rho); % C-response [m]
    % rho_a = abs(C).^2*mu_0.*omega; % Apparent resistivity [Ohm.m]
    % phi = atand(imag(C)/real(C)) + 90; % Impedance phase [deg]

    
figure
loglog(T,rho_a,'.-')
% plot(T, rho_a, '-o')
grid on
xlabel('Period T [s]','FontSize',fs)
ylabel('Apparent resistivity \rho_a [\Omega\cdotm]','FontSize',fs)
title('MT forward model with synthetic data','FontSize',fs)



