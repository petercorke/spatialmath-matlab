%H2E Homogeneous to Euclidean 
%
% E = H2E(H) is the Euclidean version (K-1xN) of the homogeneous 
% points H (KxN) where each column represents one point in P^K.
%
% See also E2H.

function e = h2e(h)

    if isvector(h)
        h = h(:);
    end
    e = h(1:end-1,:) ./ repmat(h(end,:), numrows(h)-1, 1);

