clear;clc;close all;
MYTOOLBOXROOT='E:\JJCao_Progam\visualization\JunJieCao_tool\toolbox';
addpath([MYTOOLBOXROOT '/jjcao_io/']);
addpath([MYTOOLBOXROOT '/jjcao_plot/']);
addpath([MYTOOLBOXROOT '/jjcao_interact/']);
addpath(genpath([MYTOOLBOXROOT '/jjcao_mesh/']) );

%%
%[V, F] = read_mesh('339.off');
%[V, F] = read_mesh('339noise_0.2.obj');

%[V, F] = read_mesh('Noisy.obj');

[V, F] = read_mesh('cylinder5670.1.obj');%   0.3473
%[V, F] = read_mesh('cylinder567.obj');

[Normal, NormalF] = compute_normal(V, F);
NormalF = NormalF';
fring = compute_face_ring(F);
Meas = zeros(1, 3*length(fring));
for i = 1:length(fring)
    for j = 1:length(fring{i})
        Meas(3*(i-1)+j)  = norm( NormalF(i, :) - NormalF( fring{i}(j), : ) );
    end
end

hist(Meas, 10);% 180/15这些细节被舍弃
%%
[n, xout] = hist(Meas, 10);% 180/15这些细节被舍弃
n(1:end-1) - n(2:end)
%scatter3(NormalF(1,:), NormalF(2,:), NormalF(3,:));
%%
