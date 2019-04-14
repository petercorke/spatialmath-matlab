function v = isoctave()
    v = true;

    try
        x = OCTAVE_VERSION;
    catch ME
        v = false;
    end
end
