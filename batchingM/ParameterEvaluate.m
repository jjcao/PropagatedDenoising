%%
%evaluate parameters
inputPath = 'E:\work\PropagatedDenoising\models\';
modelPath = {'Block0.4', 'Bunny0.2', 'compareJustin', ...
    'compareWang', 'compareWang', 'compareWang', ...
    'Fandisk0.3', 'Fandisk0.7', 'Julius0.2', ...
    'Nicolo0.2', 'SharpSphere0.3', 'Twelve0.5Impulsive'};
filename = {'Noisy', 'Noisy', 'trim-star_n0.5_rescaled', ...
    'bunnynoise', 'torusnoise',  'twelvenoise', ...
    'Noisy', 'Noisy','Noisy', ...
    'Noisy', 'Noisy', 'Noisy'};

% filename = 'Original';

algorithm = {'GuidedMeshNormalFiltering', 'ProjectivePropagationMeshFiltering', ...
    'ShortestPropagationMeshFiltering'};

para1 = '--faceNeighborType'; para1V = '0';
para2 = '--FaceDist'; para2V = '2';
para3 = '--SigmaS'; para3V = '1.0';
para4 = '--SigmaR'; para4V = '0.35';%para4已经没有用了
para5 = '--NormalIterNum'; para5V = '20';
para6 = '--VertexIterNum'; para6V = '10';

%% simple call
for i=1:length(filename)
    infile = [inputPath modelPath{i} '\' filename{i} '.obj'];
    for j=1:length(algorithm)        
        outdenoise = [modelPath{i} '_' filename{i} '_' algorithm{j}(1:6) '_default'];
        system(['..\src\x64\Release\Denoising.exe', ' ', infile, ' ', outdenoise, ' ', algorithm{j}]);     
    end
end

%% parameters test
j = 3;
for i=length(filename)-1:length(filename)
    infile = [inputPath modelPath{i} '\' filename{i} '.obj'];
    for k = 2:2:6 % FaceDist
        para2V = num2str(k);
        outdenoise = [modelPath{i} '_' filename{i} '_' algorithm{j}(1:6) '_faceDist_' para2V];
    % for i = 0.8:0.2:1.6 % SigmaS    
    %     para3V = num2str(i);        
    %     outdenoise = [modelPath{i} '_' filename{i} '_' algorithm{j}(1:6) '_SigmaS_' para3V];
        
        system(['..\src\x64\Release\Denoising.exe', ' ', infile, ' ', outdenoise, ' ', algorithm{j}, ' ',...
        para1, ' ', para1V, ' ', ...
        para2, ' ', para2V, ' ', ...
        para3, ' ', para3V, ' ', ...
        para4, ' ', para4V, ' ', ...
        para5, ' ', para5V, ' ', ...
        para6, ' ', para6V]); 
    end
end



