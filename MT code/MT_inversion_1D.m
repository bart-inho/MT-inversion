% Magnetotelluric (MT) 1-D inversion

clear

%% Loading of data, constants, etc...

load freq.mat % [1/s] Frequencies of measurements
load Z.mat  % [mm/s] Impedance tensor for 3 stations, with each component in:
    % - 1st dimension of array (rows) is relative to a frequency 'freq'
    % - 2nd dimension of array (columns) is relative to one component of the tensor Z
        % Z(:,1,:) := Zxx
        % Z(:,2,:) := Zxy
        % Z(:,3,:) := Zyx
        % Z(:,4,:) := Zyy
    % - 3rd dimension of array is relative to a station
load variance.mat   % [-] Variance for 3 stations, with each component in:    
    % - 1st dimension of array (rows) is relative to a frequency 'freq'
    % - 2nd dimension of array is relative to a station

T = 1./freq; % [s] Periods
ome = 2*pi.*freq; % [1/s] Angular frequency
mu0 = 4*pi*1e-7; % [kg.m.s^-2.A^-2] Magnetic permeability of free space

% Station of interest for the inversion
stn = 1;

% Impedance tensor transformed to 1D
Z_B = (Z(:,2,stn)-Z(:,3,stn))./2; % Berdichevsky average: Equation (8.8) (Simpson & Bahr, 2005)
% Z_B = (real(Z(:,2,k))-real(Z(:,3,k)))./2 + 1i* (imag(Z(:,2,k))-imag(Z(:,3,k)))./2;
Z = Z_B.*1e3; % [m/s] conversion from mm/s to m/s

% C-response
re_c = (1./ome).*imag(Z); % Eq. (???)
im_c = (-1./ome).*real(Z); % Eq. (???)
C = re_c + 1i*im_c;

rho_a = abs(C).^2*mu0.*ome; % [Ohm.m] Apparent resistivity - Eq. (2.25) from Simpson & Bahr (2005)
phi = atand(im_c./re_c)+90; % [deg] Impedance phase lag - Eq. (2.41) from Simpson & Bahr (2005)

% Creation of the structure in depth
nlayer = 21; % Number of layers

% Thicknesses of layers [m]
thick = ones(nlayer,1);
thick(1) = 50;
for j=2:nlayer-1
    thick(j) = 1.2*thick(j-1);
end, clear j
thick(end) = 60e3;

% Depths of layer interfaces [m]
z = zeros(size(thick));
for i = 1:length(thick-1)
    z(i+1) = sum(z)+thick(i);
end

% Electrical conductivity of last layer [S/m]
sigma = ones(nlayer,1); 
sigma(:) = 1/rho_a(end);


%% 1D inversion

% Tip : you need to find the best lagrange lambda parameter that is in the
% elbow of the L-curve (see Irving slides and Constable 1987). Lambda is
% different for each station.

% Derivative matrix D
sm = ones(nlayer,1);
s0 = zeros(nlayer, 1);
D = spdiags([-sm sm], -1:0, nlayer, nlayer);
D(1, :) = 0;

lamb_vec = 14;%logspace(-1, 4, 1e2)'; % Lagrange parameters
chi2_vec = zeros(size(lamb_vec)); % Chi-squared initialization
R1D = zeros(size(chi2_vec)); % Roughness parameter initialization

% Loop over the Lagrange parameters
for s = 1:length(lamb_vec)
    lambda = lamb_vec(s); % Lagrange parameter used for the inversion

    % Error matrix E
    E = diag([1./variance(:,stn); 1./variance(:,stn)]);

    % Inversion    
    thick_mod = thick(1:end-1);

    M = length(freq);
    N = length(sigma);

    m_iter = log(sigma); % Initial model (has to be log, cf. 'inversion_step.m')
    for iter = 1:25 % Number of iterations (after 25 it doesnot change a 
                    % lot anymore, but you can put more iterations if you 
                    % want to check).

        dm = 1e-4; % delta_model for numerical differentiation
        [m_iter, chi2] = inversion_step(C, T, thick_mod, m_iter, M, N, dm, E, lambda, D); % Inversion
    end
    m_end = m_iter;
    
    chi2_vec(s) = chi2;
    R1D(s) = norm(D*m_end).^2;
end

% Lagrange parameter selection
% dchi2 = abs(1-(chi2_vec-M));
% ilambda = find(dchi2==min(dchi2));
% lambda = lamb_vec(ilambda);


% Plot of results
fs = 13; % ,'FontSize',fs
fig = 21;

% Figure X1
figure(fig), clf
sgtitle(['Data station ',num2str(stn)],'FontSize',fs+2)
xLim = [min(T) max(T)];
subplot(2,1,1) % Apparent resistivity
loglog(T, rho_a)
xlabel('log_{10}T [s]')
ylabel('Apparent resistivity \rho_a [\Omega\cdotm]','FontSize',fs)
ylim([1 1e3])
xlim(xLim)
grid on
subplot(2,1,2) % Impedance phase lag
semilogx(T, phi)
xlabel('Periods T [s]','FontSize',fs)
ylabel('Phase \phi [deg]','FontSize',fs)
ylim([-180 180])
xlim(xLim)
grid on

% Figure X2
fig = fig+1;
figure(fig), clf
plot(chi2_vec-M, R1D,'+-')
hold on
% plot(chi2_vec(ilambda)-M,R1D(ilambda),'or')
title(['L-curve, station ',num2str(stn)],'FontSize',fs)
xlabel('\chi^{2}-M','FontSize',fs)
ylabel('R_{1D}','FontSize',fs)
grid on
hold off

% Figure X3
fig = fig+1;
figure(fig), clf
semilogx(1./exp(m_end), z(1:end-1)./1e3)
title(['Station ',num2str(stn)])
xlabel('Modeled resisitivies \rho [Ohm.m]')
ylabel('Depth z [km]')
grid on
axis ij

% % Figure X4
% fig = fig+1;
% figure(fig), clf
% plot(m_end,'o-')



