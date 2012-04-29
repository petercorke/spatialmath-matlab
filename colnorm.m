%COLNORM Column-wise norm of a matrix
%
% CN = COLNORM(A) is an Mx1 vector of the normals of each column of the
% matrix A which is NxM.
function n = colnorm(a)

	n = sqrt(sum(a.^2));
