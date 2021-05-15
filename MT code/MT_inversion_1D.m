%% Magnetotelluric (MT) 1-D inversion

% Section 1

clear

% Station of interest for the inversion
stn = 3;
    % (1) Station 901
    % (2) Station 902
    % (3) Station 903

disp(['Station 90',num2str(stn),'.'])
disp('Begin Section 1.')

% Loading of data, constants, etc...

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
omega = 2*pi.*freq; % [1/s] Angular frequency
mu0 = 4*pi*1e-7; % [kg.m.s^-2.A^-2] Magnetic permeability of free space

% Impedance tensor transformed to 1D
Z_B = (Z(:,2,stn)-Z(:,3,stn))./2; % Berdichevsky average: Equation (8.8) (Simpson & Bahr, 2005)
% Z_B = (real(Z(:,2,k))-real(Z(:,3,k)))./2 + 1i* (imag(Z(:,2,k))-imag(Z(:,3,k)))./2;
Z = Z_B.*1e3; % [m/s] conversion from mm/s to m/s

% C-response
re_c = (1./omega).*imag(Z); % Eq. (???)
im_c = (-1./omega).*real(Z); % Eq. (???)
C = re_c + 1i*im_c;

rho_a = abs(C).^2*mu0.*omega; % [Ohm.m] Apparent resistivity - Eq. (2.25) from Simpson & Bahr (2005)
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
z = zeros(size(thick,1)+1,1);
for i = 1:length(thick)
    z(i+1) = z(i)+thick(i);
end

% Electrical conductivity of last layer [S/m]
sigma = 1/rho_a(end) * ones(nlayer,1);


% % 1D inversion

% Tip : you need to find the best lagrange lambda parameter that is in the
% elbow of the L-curve (see Irving slides and Constable 1987). Lambda is
% different for each station.

% Derivative matrix D
sm = ones(nlayer,1);
s0 = zeros(nlayer, 1);
D = spdiags([-sm sm], -1:0, nlayer, nlayer);
D(1, :) = 0;

% Error matrix E
E = diag([1./variance(:,stn); 1./variance(:,stn)]);

thick_mod = thick(1:end-1);

% Dimensions
M = length(freq);
N = length(sigma);

% Initialization
% % % lamb_vec = logspace(-1, 2, 50)'; % Lagrange parameters 1.0235;%
lamb_vec = logspace(4, -1, 150)'; % Lagrange parameters 1.0235;%
chi2_vec = zeros(size(lamb_vec)); % Chi-squared initialization
R1D_vec = zeros(size(lamb_vec)); % Roughness parameter initialization
m_vec = zeros(length(sigma),length(lamb_vec)); % Modeled conductivities initialization

% Loop over the Lagrange parameters
for s = 1:length(lamb_vec)
    lambda = lamb_vec(s); % Lagrange parameter used for the inversion

    m_iter = log(sigma); % Initial model (has to be log, cf. 'inversion_step.m')
    for iter = 1:25 % Number of iterations (after 25 it doesnot change a 
                    % lot anymore, but you can put more iterations if you 
                    % want to check).

        dm = 1e-4; % delta_model for numerical differentiation
        [m_iter, chi2] = inversion_step(C, T, thick_mod, m_iter, M, N, dm, E, lambda, D); % Inversion
    end
    m_end = m_iter;
    
    chi2_vec(s) = chi2;
    R1D_vec(s) = norm(D*m_end).^2;
    m_vec(:,s) = m_end;
end

disp('Inversion 1D done.')


% Plot parameters
fs = 13; % ,'FontSize',fs
lw = 1.5; % ,'LineWidth',lw
fig = stn*10+1;

% L-curve
% Figure X1
figure(fig), clf
plot(chi2_vec-M, R1D_vec,'+-','LineWidth',lw)
hold on
% plot(chi2_vec(ilambda)-M,R1D(ilambda),'or')
title(['L-curve: station 90',num2str(stn)],'FontSize',fs)
xlabel('\chi^{2}-M','FontSize',fs)
ylabel('R_{1D}','FontSize',fs)
axis equal
grid on
hold off


% %%%%%%%%%%%%%%%%%%%%%%%%%%%% RETURN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('End of Section 1.')
disp('Now choose data tip (lambda) in the elbow of L-curve, then run Section 2 or set index parameter in Section 2.')
% Export cursor data to workspace from selected data tip in L-curve and
% name it 'cursor_info'.
return
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Section 2 - Lagrange parameter graphically selected

disp('Begin Section 2.')

% Plot parameters
fs = 13; % ,'FontSize',fs
lw = 1.2; % ,'LineWidth',lw

% index = cursor_info.DataIndex; % Choose data tip in elbow of L-curve
index = 0; % Set chosen lambda and skip choosing data tip step
if index==0
    if stn==1
        index = 120;
    elseif stn==2
        index = 117;
    elseif stn==3
        index = 121;
    else
        error('No station selected.')
    end
end

lambda = lamb_vec(index);
chi2 = chi2_vec(index);
m_end = m_vec(:,index);

% % % % % % % %
% Je pense que vous avez confondu la résistivité apparente rho_a qui est
% unique pour chaque mesure donc fréquence, et qui vous donne une idée sur
% la structure interne, avec les résistivités des couches rho_mod. rho_a a
% donc une dimension 36*1 alors que rho_mod 21*1
% % % % % % % %
% NEW LINE : 
rho_model = 1./exp(m_end);

disp(['Lagrange parameter lambda = ',num2str(lambda),' chosen.'])


% Forward model
% % % % % % % %
% même erreur que la dernière fois... ce qui sort de inversion_step c'est
% log(sigma) donc ill faut que vous mettiez exp dans votre input à
% Wait_recursion...
% % % % % % % %
% FAUX : [C_mod,rho_mod,phi_mod] = Wait_recursion(T,thick,1./m_end);
% NEW LINE : 
[C_forward,rho_a_forward,phi_forward] = Wait_recursion(T,thick,1./exp(m_end));

disp('Forward model done.')

% Plots

% Figure X2
figure(fig+1), clf
sgtitle(['Station 90',num2str(stn),...
    ' : \chi^{2} = ',num2str(chi2),...
    ' ; \lambda = ',num2str(round(lambda,2))],...
    'FontSize',fs+2)
padded = 0.7;
xLim = [min(T)*(1-padded) max(T)*(1+padded)];
set(gcf,'Position',[100 100 800 500])
% --- subplot 1 ---
subplot(2,2,1) % C VS T
% % % In this plot the objective if to see the results, so you use the C,
% rhoa and phi calculated with the forward pb from the restulted model
% FAUX : plot(T,re_c./1e3,'m','LineWidth',lw)
% NEW LINE :
plot(T,real(C_forward)./1e3,'m','LineWidth',lw)
hold on
% FAUX : plot(T,im_c./1e3,'g','LineWidth',lw)
% NEW LINE :
plot(T,imag(C_forward)./1e3,'g','LineWidth',lw)
xlabel('T [s]','FontSize',fs)
ylabel('C-response [km]','FontSize',fs)
xlim(xLim)
legend('real part', 'imaginary part', 'Location', 'NorthWest')
grid on
set(gca,'XScale','log');
% --- subplot 3 ---
subplot(2,2,3) % rho VS z
% % % stairs(1./exp([m_end;m_end(end)]), z./1e3,'b','LineWidth',lw)
stairs([rho_model(1); rho_model], [z]/1e3,'b','LineWidth',lw)
% title(['Model: station ',num2str(stn),' ; \lambda = ',num2str(lambda)])
xlabel('Modeled resistivity \rho [\Omega\cdotm]','FontSize',fs)
ylabel('Depth z [km]','FontSize',fs)
ylim([0 10])
set(gca,'XScale','log')
grid on
axis ij
% --- subplot 2 ---
subplot(2,2,2) % rho VS T
% % % 
loglog(T, rho_a_forward,'-b','LineWidth',lw)
hold on
% % % 
loglog(T, rho_a,'or','LineWidth',lw)
xlabel('T [s]','FontSize',fs)
ylabel('Apparent resistivity \rho_a [\Omega\cdotm]','FontSize',fs)
legend('modeled', 'observed', 'Location', 'NorthEast')
ylim([0 1e3])
xlim(xLim)
grid on
hold off
% --- subplot 4 ---
subplot(2,2,4) % phi VS T
% % % 
semilogx(T, phi_forward,'-b','LineWidth',lw)
hold on
semilogx(T, phi,'or','LineWidth',lw)
xlabel('T [s]','FontSize',fs)
ylabel('Phase \phi [deg]','FontSize',fs)
legend('modeled', 'observed', 'Location', 'SouthEast')
ylim([-180 180])
xlim(xLim)
grid on
hold off


% L-curve with selected lambda
% Figure X3 
figure(fig+2), clf
plot(chi2_vec-M, R1D_vec,'+-','LineWidth',lw)
hold on
scatter(chi2_vec(index)-M,R1D_vec(index),50,'o','LineWidth',2)
hold off
title(['L-curve: station 90',num2str(stn)],'FontSize',fs)
xlabel('\chi^{2}-M','FontSize',fs)
ylabel('R_{1D}','FontSize',fs)
axis equal
xlim([0 10])
ylim([0 10])
legend('L-curve',['\lambda = ',num2str(round(lambda,2))],'FontSize',fs)
grid on

disp('End of Section 2.')
disp('End of code.')


%% (Optional section) Saving/printing figures

disp('Begin figure printing.')

% Step 1: Make sure all figures to be saved are open

% Step 2: Set printing parameters
filedir = 'figure-inv/'; % file directory
fileformat = '-depsc'; % printing format

% Step 3: Save L-curves with chosen lambda & results of inversion
for ifig=12:10:32
    
    filestn = ['_Station90',num2str((ifig-2)/10)]; % Station number
    
    % Results of inversion
    figure(ifig)
    filename = [filedir,'Results',filestn];
    print(filename,fileformat)
    
    % L-curves
    figure(ifig+1)
    filename = [filedir,'Lcurve',filestn];
    print(filename,fileformat)
    
end

disp('Figures printed.')
disp('End.')




