function [C] = Wait_recursion(omega,z,sigma,mu_0)
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


% Initialization
q = sqrt(1i*mu_0*sigma(end)*omega); % [1/m] Inverse homogeneous half-space model transfer function 
C = 1./q; % Transfer function 

% Wait's recursion algorithm
for n=length(z):-1:2
    q = sqrt(1i*mu_0*sigma(n).*omega);
    % Wait's recursion formula
    C = (       q.*C  + tanh(q.*(z(n)-z(n-1)))  ) ./ ...
        (q.*(1+(q.*C .* tanh(q.*(z(n)-z(n-1))))));
end



