%% Create a personalized pagerank matrix from an input walk matrix.
% Params are passed in via a params.hdf5 file. Required params:
%  - Walk matrix W
%  - Page Rank dampening factor alpha
%  - Outputfile

% By default, hdf5 stores matrices in row-major order and
% MATLAB stores matrices in column-major order, so we transpose
% W' and PPR below.

W = h5read('params.hdf5','/W');
alpha = h5read('params.hdf5','/alpha');
outputfile = h5read('params.hdf5','/outputfile');

n = length(W);
PPR = (1-alpha)*inv(eye(n,n)-alpha*W)';
h5create(outputfile,'/PPR',size(PPR))
h5write(outputfile,'/PPR',PPR)
