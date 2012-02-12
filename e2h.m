%E2H Euclidean to homogeneous

function h = e2h(e)
    h = [e; ones(1,numcols(e))];
