% Dimensionality

% •	Z: les composants des tenseurs d'impédance 
%   o	colonne 1 : Zxx
%   o	colonne 2 : Zxy
%   o	colonne 3 : Zyx
%   o	colonne 4 : Zyy

% dimension 1 : frequency
% dimension 2 : impedence tensor
% dimension 3 : datasts

load Z
load freq

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

if kappa < 0.1
    if sigma < 0.05
        class(:,:) = 1; % class 1a (Simple 1-D model)
    else
        class(:,:) = 1.5; % class 1b (Simple 2-D model)
    end
else
    if mu < 0.05
        class(:,:) = 2; % class 2 (regional 1-D model with galvanic distortion)
    else
        class(:,:) = 3; % other shit
    end
end

% (5.20)
function commutator = commu(x, y)
    commutator = real(x).*imag(y)-real(y).*imag(x);
end

% S2 = Z(:,1,:);