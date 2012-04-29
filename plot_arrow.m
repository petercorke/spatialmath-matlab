%PLOT_ARROW Plot arrow
%
% PLOT_ARROW(P, OPTIONS) draws an arrow from P1 to P2 where P=[P1; P2].
%
% See also ARROW3.
function plot_arrow(p, varargin)
    mstart = p(1:end-1,:);
    mend = p(2:end,:);
    %mstart = p;
    %mend = [p(2:end,:); p(1,:)];

    arrow3(mstart, mend, varargin{:});
