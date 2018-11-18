
function tests = PluckerTest
  tests = functiontests(localfunctions);
end

function constructor_test(tc)
end

function methods_test(tc)
    % intersection
        px = Plucker([0 0 0], [1 0 0]);  % x-axis
    py = Plucker([0 0 0], [0 1 0]);  % y-axis
    px1 = Plucker([0 1 0], [1 1 0]); % offset x-axis
    
    verifyEqual(tc, px.origin_distance(), 0);
    verifyEqual(tc, px1.origin_distance(), 1);
    verifyEqual(tc, px1.origin_closesst(), [0 1 0]');

    
    
    px.intersect(px)
      px.intersect(py)  
          px.intersect(px1)
end

function intersect_test(tc)
    px = Plucker([0 0 0], [1 0 0]);  % x-axis
    py = Plucker([0 0 0], [0 1 0]);  % y-axis
    
    plane.d = [1 0 0]; plane.p = 2; % plane x=2
    
    px.intersect_plane(plane)
    py.intersect_plane(plane)
end

