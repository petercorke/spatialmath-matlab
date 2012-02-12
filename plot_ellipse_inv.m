%plot_ellipse Plot an ellipse
%
%   plot_ellipse(A, xc, ls)
%
%   ls is the standard line styles.

function h = plot_ellipse_inv(A, xc, varargin)

    if nargin == 1
        h = plot_ellipse(inv(A));
    elseif nargin == 2
        h = plot_ellipse(inv(A), xc);
    else
        h = plot_ellipse(inv(A), xc, varargin{:});
    end
