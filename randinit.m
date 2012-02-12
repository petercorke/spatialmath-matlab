function randinit(seed)

    if nargin == 0
        seed = 0;
    end

    stream = RandStream.getDefaultStream;
    stream.reset()

