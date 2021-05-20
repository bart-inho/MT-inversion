function [m_after, chi2] = inversion_step(C, T, d, m, M, N, dm, E, lambda, D)
% 
% Camp de Geophysique d'Exploration
% Projet 5: Magnetotellurique
% Bastien Ruols
% modified by Barthelemy Anhorn & Bruno Galissard de Marignac
% 
%-----------------------------------------------------------------------
% inversion_step(C, T, d, m, M, N, dm, E, lambda, D, step) gives you the
% next iteration of the model as well as the data misfit chi2.
%
% ON INPUT:
%   C : C-response from the data
%   T : vector with periods [s]
%   d : vector with layers thicknesses [m]
%   m : model as log(sigma)
%       -> sigma : electrical conductivity vector from the structure in
%                  [S.m^(-1)] (SI unit)
%   M : number of measurements by station
%   N : number of layers in the subsurface structure
%   dm : model step for numerical differentiation
%   E : error diagonal matrix of size (2*M, 2*M) to fit with Jacobian matrix
%       -> 1st half (M*M) matrix : error matrix for real part  of data
%       -> 2nd half (M*M) matrix : error matrix for imag part of data 
%   lambda : Lagrange parameter
%   D : Diagonal (N*N) matrix from Occam's method
%
% ON OUTPUT:
%   m_after : model after iteration EN LOG PTN
%   chi2 : data misfit chi2 parameter
%-----------------------------------------------------------------------

y_obs = log(C);
% C_mod = C_wait(T, d, 1./exp(m)); % Bastien Ruols
C_mod = Wait_recursion(T, d, 1./exp(m)); % Bart & Bruno
y_mod = log(C_mod);
delta_d = [real(y_obs-y_mod); imag(y_obs-y_mod)];

J = jacobian_matrix(M,N,m,dm,T,d,y_mod);

delta_m = (J'*E*J + lambda*(D'*D))\(J'*E*delta_d - lambda*(D'*D)*m);
m_after = zeros(N,1);
m_after(1:end-1) = m(1:end-1) + delta_m(1:end-1);
m_after(end) = m(end); % keep the model of the last layer fixed

chi2 = delta_d'*E*delta_d;

