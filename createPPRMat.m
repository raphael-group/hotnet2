%% Create a personalized pagerank matrix from an input walk matrix.
% Params are passed in via a params.mat file. Required params:
%  - Walk matrix W
%  - Page Rank dampening factor alpha
%  - Outputfile (must end in .mat)

load params;
I = eye( size(W) );
n = length(W);
PPR = zeros( size(W) );
Z = inv(I - alpha*W);
parfor i=1:n
    X_u = zeros(n, 1);
    X_u(i) = 1;
    PPR(i, :) = Z * ((1-alpha)*X_u);
end
save(outputfile, 'PPR');

