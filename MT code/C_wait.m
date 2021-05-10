function [C_mod, varargout] = C_wait(T, d, rho)
%                    C_mod = C_wait(T, d, rho)
% [C_mod,rhoa_mod,phi_mod] = C_wait(T, d, rho)
%
%-----------------------------------------------------------------------
% C_wait(t,rho,d) computes the complex C-response of a layered Earth
%                 using the WAIT's algorithm
%
% ON INPUT:
%   t        = vector with periods [s]
%   d(1:M-1) = vector with layer thickness in [m]
%   rho(1:M) = vector with layer resistivity in [Omega*m]
%     for a model consisting of M layers (the bottom "layer" is a uniform
%     halfspace that extends to infinity, thus its thickness is not 
%     explicitely needed.)
%
% ON OUTPUT:
%   C_mod    = complex C-response, [m]
%   rhoa_mod = apparent resistivity [Omega*m]
%   phi_mod  = phase [deg]
%-----------------------------------------------------------------------

% convert to column vectors
T = T(:);
d = d(:);
rho = rho(:); 

nlayer = length(rho);

omi    = 1i*7.895683520871486e-006./T; % i * mu_0 * 2*pi/T
C_mod  = 1./sqrt(omi/rho(nlayer)); % C-response from last layer

% for n = nlayer-1:-1:1
%     C_K   = sqrt(omi./rho(n));
%     C_mod = exp(-2*C_K*d(n)).*(C_K.*C_mod-1)./(C_K.*C_mod+1);
%     C_mod = (1+C_mod)./(1-C_mod)./C_K;
% end

for n = nlayer-1:-1:1
    q_layer   = sqrt(omi./rho(n));
    E = exp(-2*q_layer*d(n));
    B = q_layer.*C_mod;
    A = E.*(B-1)./(B+1);
    C_mod = (1+A)./(1-A)./q_layer;
end

if nargout > 1
    varargout{1} = 7.895683520871486e-6./T.*abs(C_mod).^2; % rho_a [Omega*m]
    varargout{2} = 90 + 180/pi*angle(C_mod);            % phi [degrees]
end

