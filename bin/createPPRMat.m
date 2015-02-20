%% Create a personalized pagerank matrix from an input walk matrix.
% Params are passed in via a params.hdf5 file. Required params:
%  - Walk matrix W
%  - Restart probability beta
%  - Outputfile

% By default, hdf5 stores matrices in row-major order and
% MATLAB stores matrices in column-major order, so we transpose
% W' and PPR below.

W = h5read('params.hdf5','/W');
beta = h5read('params.hdf5','/beta');
outputfile = char(h5read('params.hdf5','/outputfile'));

n = length(W);
PPR = beta*inv(eye(n,n)-(1-beta)*W)';
h5create(outputfile,'/PPR',size(PPR))
h5write(outputfile,'/PPR',PPR)
