function [C,rho_a] = Wait_recursion(omega,z,sigma,mu_0)
% Wait's recursion formula (2.33) from (Simpson & Bahr, 2005)
%
% Inputs:
% - omega: angular frequency [Hz]
% - z: depth of bottom interface [km or m ?]
% - sigma: electrical conductivity

q = sqrt(1i*mu_0.*sigma.*omega);

C = zeros(size(z));
C(end) = 1/q(end);

for n=length(z):2
    % Wait's recursion formula
    C(n-1) =     ( q(n)*C(n) + tanh(q(n)*(z(n)-z(n-1)) )) / ...
            (q(n)+(q(n)*C(n) * tanh(q(n)*(z(n)-z(n-1)) )) );
end

C_norm = sqrt(real(C).^2+imag(C).^2);

rho_a = C_norm.^2.*mu_0.*omega;

