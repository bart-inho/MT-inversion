% Magnetotelluric (MT) forward model

clear

% Parameters
load freq.mat;
freq = freq';                       % [Hz] frequency of signal
T = 1./freq;                        % [s] period of signal
omega = 2*pi.*freq;                 % [Hz] angular frequency of signal
z = [0 4e3 5e2 1e3 3e3]';           % [m] depths
                                    %     of ocean, sediments and two layers of oceanic crust
sigma = [3.2 1.5 0.03 0.004 1]';    % [S/m] electrical conductivites
                                    %     of ocean, sediments and two layers of oceanic crust

% Variables
mu_0 = 4*pi*1e-7; % [H/m] magnetic permeability of free space

% Wait's recursion formula
[C,rho_a,phi] = Wait_recursion(omega,z,sigma,mu_0);

% Regrouping data
% data = [freq';rho_a';phi'];

figure(1),clf
plot(rho_a,omega,'-o')
xlabel('Apparent resistivity \rho_a [\Omega\cdotm]')
ylabel('Angular frequency \omega [Hz]')
