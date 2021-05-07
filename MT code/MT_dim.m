% Dimensionality
% -	Z: les composants des tenseurs d'impédance 
%   -	colonne 1 : Zxx
%   -	colonne 2 : Zxy
%   -	colonne 3 : Zyx
%   -	colonne 4 : Zyy
%       - dimension 1 : Z relative to frequency
%       - dimension 2 : impedence tensor
%       - dimension 3 : datasets

load Z

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
class(kappa> 0.1 & mu>= 0.05) = 3; % other shit

figure()
heatmap(class)

% (5.20)
function commutator = commu(x, y)
    commutator = real(x).*imag(y)-real(y).*imag(x);
end

% S2 = Z(:,1,:);