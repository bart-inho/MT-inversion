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
q = zeros(size(omega)); % [m] inverse homogeneous half-space model transfer function 
q(:) = sqrt(1i*mu_0*sigma(end)*omega);

% Transfer function initialization
C = zeros(size(omega));
C(:) = 1./q(:);

% Wait's recursion algorithm
for n=length(z):-1:2
    q = sqrt(1i*mu_0*sigma(n).*omega);
    % Wait's recursion formula
    C = (      q.*C + tanh(q.*(z(n)-z(n-1)))  ) ./ ...
        (q.*(1+(q.*C .* tanh(q.*(z(n)-z(n-1))))));
end

% Apparent resistivity
rho_a = sqrt(real(C).^2+imag(C).^2)*mu_0.*omega; % [Ohm.m]
% rho_a = norm(C)^2*mu_0*omega;

% Impedance phase
phi = atand(imag(C)/real(C)) + 90; % [deg]




