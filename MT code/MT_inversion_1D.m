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

T = 1./freq; % [s] periods
ndata = length(freq);

%% set-up of the inversion

% make the data 1D
Z(:,1,:) = 0;
Z(:,4,:) = 0;
Z_B = (Z(:,2,:)-Z(:,3,:))./2; % Berdichevsky average: Equation (8.8) (Simpson & Bahr, 2005)
Z(:,2,:) =  Z_B;
Z(:,3,:) = -Z_B;

% creation of the structure in depth
nlayer = 21; % Number of layers
thick = ones(nlayer,1);
thick(1) = 50; % [m] Thickness of top layer
sigma = ones(nlayer,1);
sigma(1) = 1e-1; % [S/m] electrical conductivity of top layer
for j=2:nlayer-1
    thick(j) = 1.2*thick(j-1); % [m] Thicknesses of layers
    sigma(j) = 2*sigma(j-1); % [S/m] electrical conductivities of layers
end, clear j
thick(end+1) = 60e3; % [m] Thickness of last layer
sigma(end+1) = 2*sigma(end); % [S/m] electrical conductivity of last layer
rho = 1./sigma; % [Ohm.m] resistivities of layers
C = Wait_recursion(T,thick,rho); % C-response of data

% starting model m_0
dz = (0:50:2e4)';
nz = length(dz);
m0 = 0.02*ones(nz,1); % reference constant conductivity model (not considered if alphas=0)


%% 1D inversion

% Tip : you need to find the best lagrange lambda parameter that is in the
% elbow of the L-curve (see Irving slides and Constable 1987). Lambda is
% different for each station.

% derivative matrix D
sm = ones(ndata,1);
s0 = zeros(ndata, 1);
D = spdiags([-sm sm s0], -1:1, ndata, ndata);
D(end, :) = 0;

% smoothness
Wm = D'*D;
lambda = 1e-2;

% error matrix E
E = spdiags(1./sqrt(variance),0,size(variance, 1),size(variance, 1)); %JI: added data weighting matrix to account for noise characteristics (see inversion notes)

% inversion
m_iter = m0;
thick_mod = dz;

M = length(freq);
N = nz;

% m_iter = m_start;
for iter = 1:25 % number of iterations (after 25 it doesnot change a 
                % lot anymore, but you can put more iterations if you 
                % want to check).

    dm = 1e-4; % delta_model for numerical differentiation
    % [m_iter, chi2] = inversion_step(C, T, d_mod, m_iter, M, N, dm, E, lambda, D);
    [m_iter, chi2] = inversion_step(C, T, thick_mod, m_iter, M, N, dm, E, lambda, D);

end
m_end = m_iter;

% plot of results





