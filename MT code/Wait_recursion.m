function [C,varargout] = Wait_recursion(T,thick,rho)
% function [C] = Wait_recursion(omega,z,sigma)
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

mu0 = 4*pi*1e-7; % [H/m] magnetic permeability of free space

% Initialization
q = sqrt(1i*mu0*2*pi./(T*rho(end))); % [1/m] Inverse homogeneous half-space model transfer function 
C = 1./q; % Transfer function 

% Wait's recursion algorithm
for n=length(thick):-1:1
    q = sqrt(1i*mu0*2*pi./(T*rho(n)));
    % Wait's recursion formula
    C = (       q.*C  + tanh(q.*thick(n))  ) ./ ...
        (q.*(1+(q.*C .* tanh(q.*thick(n)))));
end

if nargout > 1
    varargout{1} = abs(C).^2*mu0.*2*pi./T; % Apparent resistivity [Ohm.m]
    varargout{2} = atand(imag(C)./real(C)) + 90; % Impedance phase [deg]
end
