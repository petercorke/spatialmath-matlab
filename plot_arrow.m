function plot_arrow(p, varargin)
    mstart = p(1:end-1,:);
    mend = p(2:end,:);
    %mstart = p;
    %mend = [p(2:end,:); p(1,:)];

    arrow3(mstart, mend, varargin{:});
