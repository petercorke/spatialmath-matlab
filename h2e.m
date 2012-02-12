%H2E Homogeneous to Euclidean 

function e = h2e(h)

    if isvector(h)
        h = h(:);
    end
    e = h(1:end-1,:) ./ repmat(h(end,:), numrows(h)-1, 1);

