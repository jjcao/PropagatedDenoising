%%
%evaluate parameters

%noisefile = 'noise339_0.2.obj';
noisefile = 'torusnoise.obj';

% algorithm = 'ProjectivePropagationMeshFiltering';
% para1 = '--faceNeighborType'; para1V = '3';

algorithm = 'ShortestPropagationMeshFiltering';
para1 = '--faceNeighborType'; para1V = '2';


para2 = '--FaceDist'; para2V = '3';
para3 = '--SigmaS';%para3V = '0.1';
para4 = '--SigmaR';%para4V = '0.5';%para4已经没有用了
para5 = '--NormalIterNum'; para5V = '30';
para6 = '--VertexIterNum'; para6V = '15';

para4V = '0.5';
%projective
%for i = 0.01:0.05:1 
%shortest
for i = 0.0001:0.005:0.2
    para3V = num2str(i);
        outdenoise = ['denoise_' num2str(i)];%how to ouput many?
        system(['Denoising.exe', ' ', noisefile, ' ', outdenoise, ' ', algorithm, ' ',...
            para1, ' ', para1V, ' ', ...
            para2, ' ', para2V, ' ', ...
            para3, ' ', para3V, ' ', ...
            para4, ' ', para4V, ' ', ...
            para5, ' ', para5V, ' ', ...
            para6, ' ', para6V]);   
end


