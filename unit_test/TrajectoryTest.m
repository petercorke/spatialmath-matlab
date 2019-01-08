%% This is for testing the Trajectory Generation functions in the robotics Toolbox
function tests = TrajectoryTest
  tests = functiontests(localfunctions);
  clc
end

function teardownOnce(tc)
    close all
end

function trinterp2_test(tc)
    T0 = transl2(1, 2) * trot2(20, 'deg');
    T1 = transl2(3, 6) * trot2(60, 'deg');
    Tm = transl2(2, 4) * trot2(40, 'deg');
    
    T = trinterp2(T0, T1, [0.5 0 1]);
    tc.verifyEqual(size(T), [3 3 3]);
    tc.verifyEqual(T(:,:,1), Tm, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,2), T0, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,3), T1, 'abstol', 1e-10);
end

%%TODO trinterp

