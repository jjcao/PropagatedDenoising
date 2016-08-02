%%
%evaluate parameters
inputPath = 'E:\work\PropagatedDenoising\models\';
modelPath = 'Block0.4';
modelPath = 'Bunny0.2';
% modelPath = 'compareWang';
modelPath = 'compareJustin';
% modelPath = 'Fandisk0.3';
modelPath = 'Fandisk0.7';
% modelPath = 'Nicolo0.2';
% modelPath = 'SharpSphere0.3';
% modelPath = 'Twelve0.5Impulsive';


filename = 'Noisy';
% filename = 'Original';
% filename = 'torusnoise';
% filename = 'trim-star_n0.5_rescaled';
% filename = 'twelvenoise';

infile = [inputPath modelPath '\' filename '.obj'];

% algorithm = 'ProjectivePropagationMeshFiltering';
% para1 = '--faceNeighborType'; para1V = '3';

algorithm = 'ShortestPropagationMeshFiltering';
para1 = '--faceNeighborType'; para1V = '0';

para2 = '--FaceDist'; para2V = '4';
para3 = '--SigmaS'; para3V = '1.2';
para4 = '--SigmaR'; para4V = '0.5';%para4已经没有用了
para5 = '--NormalIterNum'; para5V = '30';
para6 = '--VertexIterNum'; para6V = '15';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%
%projective
%for i = 0.01:0.05:1 
%     para3V = num2str(i);
%shortest
for i = 3:2:8 % FaceDist
    para2V = num2str(i);
% for i = 0.8:0.2:1.6 % SigmaS    
%     para3V = num2str(i);
    outdenoise = [modelPath '_' filename '_faceDist_' num2str(i)];%how to ouput many?
    system(['..\src\x64\Release\Denoising.exe', ' ', infile, ' ', outdenoise, ' ', algorithm, ' ',...
        para1, ' ', para1V, ' ', ...
        para2, ' ', para2V, ' ', ...
        para3, ' ', para3V, ' ', ...
        para4, ' ', para4V, ' ', ...
        para5, ' ', para5V, ' ', ...
        para6, ' ', para6V]);   
end


