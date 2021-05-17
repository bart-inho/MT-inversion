% Magnetotelluric (MT) and Dimensionality
%
% Sources: (Simpson & Bahr, 2005)

clear

load freq.mat % [1/s] Frequencies of measurements
load Z.mat  % [mm/s] Impedance tensor for 3 stations, with each component in:
    % - 1st dimension of array (rows) is relative to a frequency 'freq'
    % - 2nd dimension of array (columns) is relative to one component of the tensor Z
        % Z(:,1,:) := Zxx
        % Z(:,2,:) := Zxy
        % Z(:,3,:) := Zyx
        % Z(:,4,:) := Zyy
    % - 3rd dimension of array is relative to a station

% -------------------------
% Verification with 1D tensor
% Impedance tensor transformed to 1D
% Z_B = (Z(:,2,:)-Z(:,3,:))./2; % Berdichevsky average: Equation (8.8) (Simpson & Bahr, 2005)
% Z(:,1,:) = 0;
% Z(:,4,:) = 0;
% Z(:,2,:) = Z_B;
% Z(:,3,:) = -Z_B;
% -------------------------

% equation (5.6)
S1 = Z(:,1,:) + Z(:,4,:);
S2 = Z(:,2,:) + Z(:,3,:);
D1 = Z(:,1,:) - Z(:,4,:);
D2 = Z(:,2,:) - Z(:,3,:);

% equation (5.10)
sigma = (D1.^2+S2.^2)./D2.^2;

% equation (5.16)
kappa = abs(S1)./abs(D2);

% equation (5.21)
mu = sqrt(commu(D1, S2) + commu(S1, D2))./abs(D2);

class = zeros(size(Z,1), size(Z,3));        

class(kappa< 0.1 & sigma< 0.05) = 1; % class 1a (Simple 1-D model)
class(kappa< 0.1 & sigma>= 0.05) = 1.5; % class 1b (Simple 2-D model)
class(kappa> 0.1 & mu < 0.05) = 2; % class 2 (regional 1-D model with galvanic distortion)
class(kappa> 0.1 & mu>= 0.05) = 3; % other

dataset = {'Station 901','Station 902','Station 903'};
classname = {'Class 1a','Class 1b','Class 2','Other complicated classes'};
% classname = cell(size(class));
% classname(class==1  ) = {'1a'};
% classname(class==1.5) = {'1b'};
% classname(class==2  ) = {'2 '};
% classname(class==3  ) = {'other'};

fs = 10; % ,'FontSize',fs

figure
h = heatmap(dataset,round(1./freq,2,'significant'),class,'FontSize',fs);
ylabel('Periods [s]')
colorbar off

% (5.20)
function commutator = commu(x, y)
    commutator = real(x).*imag(y)-real(y).*imag(x);
end

% S2 = Z(:,1,:);