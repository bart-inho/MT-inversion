% Magnetotelluric (MT)
% Homogeneous half-space
% 
% Camp de Geophysique d'Exploration
% Projet 5: Magnetotellurique
% Barthelemy Anhorn & Bruno Galissard de Marignac
% 
% Second exercise


% Rows: relative to medium
%       ocean, crust, upper mantle, lower mantle, core
% Columns: relative to frequence

clear

fs = 13;
mu_0 = 4*pi*1e-7; % [H/m] magnetic permeability of free space

T = linspace(1,1e4*3600*24*365.25); % [s] period
omega = 2*pi./T; % [1/s] angular frequency
sigma = [3.2;1e-2;5e-3;5;5e5]; % [S/m] electrical conductivities
q = 1i*mu_0*sigma*omega;
skin = 1./(real(sqrt(q))); % [m] skin-depth

E1 = 1; % 
t = linspace(0,1e6,length(omega));

E_0 = E_field(E1,omega,t,q,0);
E_skin = E_field(E1,omega,t,q,skin);


figure
plot(T,skin*1e-3,'.-')
grid on
xlabel('Period T [s]','FontSize',fs)
ylabel('Skin depth [km]','FontSize',fs)
title('MT homogeneous half-space','FontSize',fs)
legend('ocean','crust','upper mantle','lower mantle','core','Location','best')

figure
plot(E_0-E_skin,'.-')
xlabel('Re(E_0-E_{skin})','FontSize',fs)
ylabel('Im(E_0-E_{skin})','FontSize',fs)
grid on


function [E] = E_field(E1,omega,t,q,z)
    E = E1*exp(1i*omega.*t-q.*z);
end

