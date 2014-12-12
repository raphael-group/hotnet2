%% Create a personalized pagerank matrix from an input walk matrix.
% Params are passed in via a params.mat file. Required params:
%  - Walk matrix W
%  - Page Rank dampening factor alpha
%  - Outputfile (must end in .mat)

load params;
n = length(W);
PPR = (1-alpha)*inv(eye(n,n)-alpha*W')
save(outputfile, 'PPR');

