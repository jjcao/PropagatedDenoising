%% 
%
% Copyright (c) 2016 Junjie Cao
clear;clc;close all;
MYTOOLBOXROOT='E:\jjcao_code_git\toolbox';
addpath([MYTOOLBOXROOT '/geom3d/geom3d']);

%% compute distance of two normal vector with angle = pi/? (in radians)
angle = pi/40:pi/40:pi/2;
sigma = zeros(size(angle));
gauss = zeros(size(angle));
n1 = [1 0 0];
for i=1:length(angle)    
    R = createRotationOy(angle(i));
    n2 = transformPoint3d(n1, R);
    dist = norm(n1-n2);
    sigma(i) = 0.25*dist;
    gauss(i) = exp(-0.5*dist^2/sigma(i)^2);
end
%%