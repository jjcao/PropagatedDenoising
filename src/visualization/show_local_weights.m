%% 
%
% Copyright (c) 2016 Junjie Cao
clear;clc;close all;
MYTOOLBOXROOT='E:\jjcao_code_git\toolbox';
addpath([MYTOOLBOXROOT '/jjcao_io/']);
addpath([MYTOOLBOXROOT '/jjcao_plot/']);
addpath([MYTOOLBOXROOT '/jjcao_interact/']);
addpath(genpath([MYTOOLBOXROOT '/jjcao_mesh/']) );

% mex -largeArrayDims ../Algorithms/dijkstra_mex.cpp ../Algorithms/dijkstra.cpp

% [M.verts, M.faces] = read_mesh('E:\work\PropagatedDenoising\models\compareWang/torusnoise.obj');
[M.verts, M.faces] = read_mesh('E:\work\PropagatedDenoising\models\Fandisk0.7/Original.obj');

[normal,normalf] = compute_normal(M.verts, M.faces);
normalf = normalf';
[A,faceCenter] = compute_dual_graph(M.faces,M.verts);
% figure('Name','Dual Mesh'); set(gcf,'color','white');hold on;
% options.ps=10;
% h = plot_graph(A, faceCenter, options);
% axis off;    axis equal; 
% camorbit(0,0,'camera'); axis vis3d; view(-90, 0);view3d rot;
faceCenter = faceCenter';
faceColor = zeros(size(M.faces));
faceColor(:,2) = 1;

%% compute and show k-ring face nneighbors
faceIdx = 1; 
faceRing = 5;
B = A;
faceNeighbors = find(B(faceIdx, :)>0);
for i = 1:(faceRing-1)
    B = B * A;
    B(B>0) = 1;
    faceNeighbors = union(faceNeighbors, find(B(faceIdx, :)>0));
end
faceNeighbors = union(faceIdx, faceNeighbors);


faceColor(faceNeighbors, 1) = 1;  
faceColor(faceIdx, :) = [1 0 0]; 

figure('Name','Input'); set(gcf,'color','white');
% h=trisurf(M.faces, M.verts(:,1), M.verts(:,2), M.verts(:,3), 'FaceColor', 'cyan', 'edgecolor',[1,0,0]);
h=trisurf(M.faces, M.verts(:,1), M.verts(:,2), M.verts(:,3), 'FaceVertexCData', faceColor);%,'edgecolor','none'
% delete(h);
axis off;    axis equal;   set(gcf,'Renderer','OpenGL'); 
camorbit(0,0,'camera'); axis vis3d; view(-90, 0); mouse3d; hold on;
% shading faceted; %shading interp;

%% compute local graph of k-ring face nneighbors
fring = compute_face_ring(M.faces);
localA = zeros(length(faceNeighbors));
for i=1:length(faceNeighbors)
    neiIdx = fring{faceNeighbors(i)}; % neighbor index
    for j = 1:length(neiIdx)
        tmp = norm(normalf(neiIdx(j),:) - normalf(faceNeighbors(i),:));
%         tmp = norm(faceCenter(neiIdx(j),:) - faceCenter(faceNeighbors(i),:));
        w = max(tmp, 0.00001);
            
        localNeiIdx = find(faceNeighbors==neiIdx(j));
        localA(i, localNeiIdx) = w;
    end
end

%% dijkstra path
localFaceIdx = find(faceNeighbors==faceIdx);
% [D,S] = perform_dijkstra(localA, localFaceIdx);
D = dijkstra_mex(localA, localFaceIdx);
localB = localA;
localB(localB>0)=1;
options.end_points = find(sum(localB,2) <3 );
localPath = perform_dijkstra_path_extraction(localA,D,options.end_points); 

faceColor = zeros(size(M.faces,1),1)+2*max(D);
faceColor(faceNeighbors) = D;
figure('Name','shortest path by Dijkstra'); set(gcf,'color','white');
h=trisurf(M.faces, M.verts(:,1), M.verts(:,2), M.verts(:,3), 'FaceVertexCData', faceColor, 'faceAlpha', 0.8);%, 'edgecolor','none'
% h=trisurf(M.faces(faceNeighbors,:), M.verts(:,1), M.verts(:,2), M.verts(:,3), 'FaceVertexCData', D);%, 'edgecolor','none'
axis off;    axis equal;  mouse3d; hold on; colorbar;
% set(h, 'faceAlpha', 0.8);

colors = [1 0 0; 0 1 0;];
for i = 1:length(localPath)
    p = faceNeighbors(localPath{i});
    plot3(faceCenter(p,1),faceCenter(p,2),faceCenter(p,3), 'color', colors(1,:));
    scatter3(faceCenter(p,1), faceCenter(p,2),faceCenter(p,3),40,'MarkerFaceColor',colors(1,:))
end

%%
p = faceNeighbors(localPath{1});
faceColor = zeros(size(M.faces,1),1)+2*max(D);
faceColor(faceNeighbors) = D;
faceColor(p) = max(D)*3;
faceColor(faceIdx) = max(D)*4;
figure('Name','shortest path'); set(gcf,'color','white');
h=trisurf(M.faces, M.verts(:,1), M.verts(:,2), M.verts(:,3), 'FaceVertexCData', faceColor, 'faceAlpha', 0.8);%, 'edgecolor','none'
% h=trisurf(M.faces(faceNeighbors,:), M.verts(:,1), M.verts(:,2), M.verts(:,3), 'FaceVertexCData', D);%, 'edgecolor','none'
axis off;    axis equal; mouse3d; hold on;
%% weight


