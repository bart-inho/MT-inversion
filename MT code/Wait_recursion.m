function [C,rho_a,phi] = Wait_recursion(omega,z,sigma,mu_0)
% Wait's recursion formula (2.33) from (Simpson & Bahr, 2005)
%
% Inputs:
% - omega: angular frequency [Hz]
% - z: depth of bottom interface [m]
% - sigma: electrical conductivity [S/m]
%
% Outputs:
% - C: Schmucker-Weidelt transfer function (Weidelt, 1972; Schmucker, 1973) [1/m]
% - rho_a: apparent resistivity [Ohm.m]
% - phi: impedance phase [deg]

% Inverse transfer function initialization
q = zeros(size(z,1),size(omega,2)); % [m] inverse homogeneous half-space model transfer function 
q(end,:) = sqrt(1i*mu_0.*sigma(end).*omega);

% Transfer function initialization
C = zeros(size(q));
C(end,:) = 1./q(end,:);

% Wait's recursion algorithm
for n=length(z):2
    q(n,:) = sqrt(1i*mu_0*sigma(n)*omega);
    % Wait's recursion formula
    C(n-1,:) = (              q(n,:).*C(n,:) + tanh( q(n,:)*(z(n)-z(n-1))    )) / ...
               ( q(n,:) * (1+(q(n,:).*C(n,:) * tanh( q(n,:)*(z(n)-z(n-1)) )) ));
end

% Apparent resistivity
% rho_a = (real(C).^2+imag(C).^2)*mu_0*omega; % [Ohm.m]
rho_a = norm(C)^2*mu_0*omega;

% Impedance phase
phi = atand(imag(C)./real(C)) + 90; % [deg]
