%E2H Euclidean to homogeneous
%
% H = E2H(E) is the homogeneous version (K+1xN) of the Euclidean 
% points E (KxN) where each column represents one point in R^K.
%
% See also H2E.

function h = e2h(e)
    h = [e; ones(1,numcols(e))];
