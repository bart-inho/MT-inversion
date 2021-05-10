% Magnetotelluric (MT) 1-D inversion

clear

%% loading of data, constants, etc...

load variance.mat
load freq.mat
load Z.mat
% Z(:,1,:) := Zxx
% Z(:,2,:) := Zxy
% Z(:,3,:) := Zyx
% Z(:,4,:) := Zyy

%% set-up of the inversion

% make the data 1D
Z(:,1,:) = 0;
Z(:,4,:) = 0;
Z_B = (Z(:,2,:)-Z(:,3,:))./2; % Berdichevsky average Equation (8.8) (Simpson & Bahr, 2005)
Z(:,2,:) =  Z_B;
Z(:,3,:) = -Z_B;

% creation of the structure in depth

% starting model m_0

%% 1D inversion

% Tip : you need to find the best lagrange lambda parameter that is in the
% elbow of the L-curve (see Irving slides and Constable 1987). Lambda is
% different for each station.

% derivative matrix D

% error matrix E

% inversion
    m_iter = m_start;
    for iter = 1:25 % number of iterations (after 25 it doesnot change a 
                    % lot anymore, but you can put more iterations if you 
                    % want to check).
        
        dm = 1e-4; % delta_model for numerical differentiation
        [m_iter, chi2] = inversion_step(C, T, d_mod, m_iter, M, N, dm, E, lambda, D);
        
    end
    m_end = m_iter;
    
% plot of results