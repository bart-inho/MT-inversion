function [J] = jacobian_matrix(M,N,m,dm,T,d,y_m)

%-----------------------------------------------------------------------
% jacobian_matrix(N,m) computes the jacobian matrix needed for Occam's 1D
%                      inversion of MT data
%
% ON INPUT:
%   M : number of measurements by station
%   N : number of layers in the subsurface structure
%   m : model as log(sigma) 
%       -> sigma : electrical conductivity vector from the structure in
%                  [S.m^(-1)] (SI unit)
%   dm : model step for numerical differentiation
%   T : vector with periods [s]
%   d : vector with layer thickness in [m]
%   y_m : log(initial complex C_response)
%
% ON OUTPUT:
%   J : jacobian matrix of length (size (2*M,N))
%       -> lines 1 to M : real parts
%       -> lines M+1 to 2*M : imaginary parts
%-----------------------------------------------------------------------

J = zeros(2*M, N);

for k = 1:N
    m_int = m;
    m_int(k) = m_int(k)+dm;
    C_mod_int = C_wait(T, d, 1./exp(m_int));
    y_m_int = log(C_mod_int);
    J(:,k) = [real(y_m_int-y_m)./dm; imag(y_m_int-y_m)./dm];
end