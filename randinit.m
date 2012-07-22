%RANDINIT Reset random number generator
%
% RANDINIT reset the defaul random number stream.
%
% See also RandStream.

function randinit(seed)

    stream = RandStream.getGlobalStream;
    stream.reset()

