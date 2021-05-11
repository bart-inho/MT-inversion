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
% unité Z : mm/s

T = 1./freq; % [s] periods
ndata = length(freq);

% set-up of the inversion
k = 1;
% make the data 1D
Z(:,1,k) = 0;
Z(:,4,k) = 0;
Z_B = (Z(:,2,k)-Z(:,3,k))./2; % Berdichevsky average: Equation (8.8) (Simpson & Bahr, 2005)
Z(:,2,k) =  Z_B;
Z(:,3,k) = -Z_B;
Z = Z_B.*1e3;

% C response
ome = 2*pi.*freq;
re_c = (1./ome).*imag(Z);
im_c = (-1./ome).*real(Z);
C_base = re_c + 1i*im_c;

rho_a = abs(C_base).^2*4*pi*1e-7.*ome;

% faire un joli plot
figure(1), clf
semilogy(log10(T), rho_a)

phase = atand(im_c./re_c)+90;
figure(2), clf
plot(log10(1./freq), phase)

% creation of the structure in depth
nlayer = 21; % Number of layers
thick = ones(nlayer,1);
thick(1) = 50; % [m] Thickness of top layer
sigma = ones(nlayer,1);
% sigma(1) = 1e-1; % [S/m] electrical conductivity of top layer
for j=2:nlayer-1
    thick(j) = 1.2*thick(j-1); % [m] Thicknesses of layers
%     sigma(j) = 2*sigma(j-1); % [S/m] electrical conductivities of layers
end, clear j
thick(end) = 60e3; % [m] Thickness of last layer
sigma(:) = 1/rho_a(end); % [S/m] electrical conductivity of last layer
% rho = 1./sigma; % [Ohm.m] resistivities of layers
% C = Wait_recursion(T,thick,rho); % C-response of data

% starting model m_0
dz = (0:50:2e4)';
nz = length(dz);
m0 = 0.02*ones(nz,1); % reference constant conductivity model (not considered if alphas=0)
%% 1D inversion

% Tip : you need to find the best lagrange lambda parameter that is in the
% elbow of the L-curve (see Irving slides and Constable 1987). Lambda is
% different for each station.

% derivative matrix D
sm = ones(nlayer,1);
s0 = zeros(nlayer, 1);

D = spdiags([-sm sm], -1:0, nlayer, nlayer);
D(1, :) = 0;
% full(D)

% smoothness
% Wm = D'*D;
lamb_vec = logspace(-1, 4, 1e2)';

chi2_tot = [];
R1d = [];

for s = 1:length(lamb_vec)
    lambda = lamb_vec(s);
    % lambda = 1;

    % error matrix E
    % E = spdiags(1./sqrt(variance),0,size(variance, 1),size(variance, 1));
    % variance -> pas std
    E = diag([1./variance(:, 1); 1./variance(:, 1)]);

    % inversion
    m_iter = sigma;
    thick_mod = thick(1:end-1);

    M = length(freq);
    N = length(sigma);

    % m_iter = m_start;

    for iter = 1:25 % number of iterations (after 25 it doesnot change a 
                    % lot anymore, but you can put more iterations if you 
                    % want to check).

        dm = 1e-4; % delta_model for numerical differentiation
        % [m_iter, chi2] = inversion_step(C, T, d_mod, m_iter, M, N, dm, E, lambda, D);
        [m_iter, chi2] = inversion_step(C_base, T, thick_mod, m_iter, M, N, dm, E, lambda, D);

    end
    m_end = m_iter;
    
    chi2_tot = [chi2_tot;chi2];
    R1d = [R1d;norm(D*m_end).^2];
end

% plot of results
z = zeros(length(thick), 1);
for i = 1:length(thick)
    z(i) = sum(thick(1:i));
end

figure(3), clf
plot(chi2_tot-M, R1d)
    
% figure(4), clf
% plot(m_end, z)
