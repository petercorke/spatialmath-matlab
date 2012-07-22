%MLABEL	labels for mplot style graph
%
%	mlabel({lab1 lab2 lab3})

function mlabel(lab, varargin)

	% find all child axes (subplots)
	h = findobj(gcf, 'Type', 'axes');

	for i=1:length(h),

		if strcmp( get(h(i), 'visible'), 'on'),
			axes(h(i))
			% get subplot number from user data (I don't know who
			% sets this but its very useful)
			sp = get(h(i), 'UserData');
			if sp == 1,
				topplot = sp;
			end
			ylabel(lab{sp}, varargin{:});
		end
	end

    if 0
        if nargin > 1,
            axes(h(topplot));	% top plot
            title(tit);
        end
    end
