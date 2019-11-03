%% This is for testing the Homogeneous Transformation functions in the robotics Toolbox

function tests = TransformationsTest
    tests = functiontests(localfunctions);
    
    clc
end

function teardownOnce(tc)
    close all
end

%% first of all check we can tell a good matrix from a bad one
function isrot_test(tc)
    R1 = diag([1 1 1]);    % proper
    R2 = diag([1 1 -1]);   % not proper
    R3 = diag([1 2 1]);    % not proper
    R4 = diag([2 0.5 1]);  % not proper
    
    % test shapes
    tc.verifyFalse( isrot(1) )
    tc.verifyFalse( isrot( zeros(2,2) ) )
    tc.verifyFalse( isrot( zeros(4,4) ) )
    tc.verifyFalse( isrot( zeros(3,1) ) )
    tc.verifyFalse( isrot( zeros(1,3) ) )
    
    % test shapes with validity check
    tc.verifyFalse( isrot(1, 1) )
    tc.verifyFalse( isrot( zeros(2,2), 1 ) )
    tc.verifyFalse( isrot( zeros(4,4) ), 1 )
    tc.verifyFalse( isrot( zeros(4,1) ), 1 )
    tc.verifyFalse( isrot( zeros(1,4) ), 1 )
    
    % test 3x3
    tc.verifyTrue( isrot(R1) )
    tc.verifyTrue( isrot(R2) )
    tc.verifyTrue( isrot(R3) )
    tc.verifyTrue( isrot(R4) )
    
    % test 3x3 with validity check
    tc.verifyTrue( isrot(R1, 1) )
    tc.verifyFalse( isrot(R2, 1) )
    tc.verifyFalse( isrot(R3, 1) )
    tc.verifyFalse( isrot(R4, 1) )
    
    % vector case
    tc.verifyTrue( isrot(cat(3, R1, R1, R1)) )
    tc.verifyTrue( isrot(cat(3, R1, R1, R1), 1) )
    tc.verifyTrue( isrot(cat(3, R1, R2, R3)) )
    tc.verifyFalse( isrot(cat(3, R1, R2, R3), 1) )
end

function ishomog_test(tc)
    T1 = diag([1 1 1 1]);    % proper
    T2 = diag([1 1 -1 1]);   % not proper
    T3 = diag([1 2 1 1]);    % not proper
    T4 = diag([2 0.5 1 1]);  % not proper
    T5 = diag([1 1 1 0]);    % not proper
    
    
    % test shapes
    tc.verifyFalse( ishomog(1) )
    tc.verifyFalse( ishomog( zeros(2,2) ) )
    tc.verifyFalse( ishomog( zeros(3,3) ) )
    tc.verifyFalse( ishomog( zeros(4,1) ) )
    tc.verifyFalse( ishomog( zeros(1,4) ) )
    
    % test shapes with validity check
    tc.verifyFalse( ishomog(1, 1) )
    tc.verifyFalse( ishomog( zeros(2,2), 1 ) )
    tc.verifyFalse( ishomog( zeros(3,3) ), 1 )
    tc.verifyFalse( ishomog( zeros(4,1) ), 1 )
    tc.verifyFalse( ishomog( zeros(1,4) ), 1 )
    
    % test 4x4
    tc.verifyTrue( ishomog(T1) )
    tc.verifyTrue( ishomog(T2) )
    tc.verifyTrue( ishomog(T3) )
    tc.verifyTrue( ishomog(T4) )
    tc.verifyTrue( ishomog(T5) )
    
    
    % test 4x4 with validity check
    tc.verifyTrue( ishomog(T1, 1) )
    tc.verifyFalse( ishomog(T2, 1) )
    tc.verifyFalse( ishomog(T3, 1) )
    tc.verifyFalse( ishomog(T4, 1) )
    tc.verifyFalse( ishomog(T5, 1) )
    
    
    % vector case
    tc.verifyTrue( ishomog(cat(3, T1, T1, T1)) )
    tc.verifyTrue( ishomog(cat(3, T1, T1, T1), 1) )
    tc.verifyTrue( ishomog(cat(3, T1, T2, T3)) )
    tc.verifyFalse( ishomog(cat(3, T1, T2, T3), 1) )
end


%% can we convert between rotation matrices and homogeneous coordinate matrices

function r2t_test(tc)
    
    % SO(3) case
    R = [1 2 3;4 5 6; 7 8 9];
    tc.verifyEqual(r2t(R),...
        [1 2 3 0; 4 5 6 0; 7 8 9 0; 0 0 0 1],'absTol',1e-10);
    
    % sequence case
    Rs = cat(3, R, R, R, R, R);
    Ts = r2t(Rs);
    verifySize(tc, Ts, [4 4 5]);
    tc.verifyEqual(Ts(:,:,2), ...
        [1 2 3 0; 4 5 6 0; 7 8 9 0; 0 0 0 1],'absTol',1e-10);
    
end


function t2r_test(tc)
    %Unit test for r2t with variables eul2tr([.1, .2, .3])
    
    % SE(3) case
    T = [1 2 3 4; 5 6 7 8; 9 10 11 12; 0 0 0 1];
    tc.verifyEqual(t2r(T),...
        [1 2 3; 5 6 7; 9 10 11],'absTol',1e-10);
    
    % sequence case
    Ts = cat(3, T, T, T, T, T);
    Rs = t2r(Ts);
    verifySize(tc, Rs, [3 3 5]);
    tc.verifyEqual(Rs(:,:,2), ...
        [1 2 3; 5 6 7; 9 10 11],'absTol',1e-10);
    
end

function rt2tr_test(tc)
    
    R = [1 2 3;4 5 6; 7 8 9];
    t = [-10; -11; -12];
    
    tc.verifyEqual(rt2tr(R, t),...
        [1 2 3 -10; 4 5 6 -11; 7 8 9 -12; 0 0 0 1],'absTol',1e-10);
    
    % sequence case
    Rs = cat(3, R, 2*R, 3*R);
    ts = cat(2, t, 2*t, 3*t);
    Ts = rt2tr(Rs, ts);
    verifySize(tc, Ts, [4 4 3]);
    tc.verifyEqual(Ts(:,:,1), ...
        [1 2 3 -10; 4 5 6 -11; 7 8 9 -12; 0 0 0 1],'absTol',1e-10);
    tc.verifyEqual(Ts(:,:,2), ...
        [2*[1 2 3 -10; 4 5 6 -11; 7 8 9 -12]; 0 0 0 1],'absTol',1e-10);
    tc.verifyEqual(Ts(:,:,3), ...
        [3*[1 2 3 -10; 4 5 6 -11; 7 8 9 -12]; 0 0 0 1],'absTol',1e-10);
end


function tr2rt_test(tc)
    %Unit test for r2t with variables eul2tr([.1, .2, .3])
    
    %% SE(3) case
    T = [1 2 3 4; 5 6 7 8; 9 10 11 12; 0 0 0 1];
    [R,t] = tr2rt(T);
    tc.verifyEqual(R, [1 2 3; 5 6 7; 9 10 11], 'absTol',1e-10);
    tc.verifyEqual(t, [4; 8; 12], 'absTol',1e-10);
    
    % sequence case
    Ts = cat(3, T, T, T, T, T);
    [Rs,ts] = tr2rt(Ts);
    verifySize(tc, Rs, [3 3 5]);
    verifySize(tc, ts, [5 3]);
    
    tc.verifyEqual(Rs(:,:,2), [1 2 3; 5 6 7; 9 10 11], 'absTol',1e-10);
    tc.verifyEqual(ts(2,:), [4 8 12], 'absTol',1e-10);
    
end

% Primitives
function rotx_test(tc)
    tc.verifyEqual(rotx(0), eye(3,3),'absTol',1e-10);
    tc.verifyEqual(rotx(pi/2), [1 0 0; 0 0 -1; 0 1 0],'absTol',1e-10);
    tc.verifyEqual(rotx(pi), [1 0 0; 0 -1 0; 0 0 -1],'absTol',1e-10);
    
    tc.verifyEqual(rotx(90, 'deg'), [1 0 0; 0 0 -1; 0 1 0],'absTol',1e-10);
    tc.verifyEqual(rotx(180, 'deg'), [1 0 0; 0 -1 0; 0 0 -1],'absTol',1e-10);
    
    syms q
    R = rotx(q);
    verifyInstanceOf(tc, R, 'sym');
    verifySize(tc, R, [3 3]);
    tc.verifyEqual(simplify(det(R)), sym(1));
    
    %test for non-scalar input
    verifyError(tc, @()rotx([1 2 3]),'SMTB:rotx:badarg');
end

function roty_test(tc)
    tc.verifyEqual(roty(0), eye(3,3),'absTol',1e-10);
    tc.verifyEqual(roty(pi/2), [0 0 1; 0 1 0; -1 0 0],'absTol',1e-10);
    tc.verifyEqual(roty(pi), [-1 0 0; 0 1 0; 0 0 -1],'absTol',1e-10);
    
    tc.verifyEqual(roty(90, 'deg'), [0 0 1; 0 1 0; -1 0 0],'absTol',1e-10);
    tc.verifyEqual(roty(180, 'deg'), [-1 0 0; 0 1 0; 0 0 -1],'absTol',1e-10);
    
    syms q
    R = roty(q);
    verifyInstanceOf(tc, R, 'sym');
    verifySize(tc, R, [3 3]);
    tc.verifyEqual(simplify(det(R)), sym(1));
    
    %test for non-scalar input
    verifyError(tc, @()roty([1 2 3]),'SMTB:roty:badarg');
end

function rotz_test(tc)
    tc.verifyEqual(rotz(0), eye(3,3),'absTol',1e-10);
    tc.verifyEqual(rotz(pi/2), [0 -1 0; 1 0 0; 0 0 1],'absTol',1e-10);
    tc.verifyEqual(rotz(pi), [-1 0 0; 0 -1 0; 0 0 1],'absTol',1e-10);
    
    tc.verifyEqual(rotz(90, 'deg'), [0 -1 0; 1 0 0; 0 0 1],'absTol',1e-10);
    tc.verifyEqual(rotz(180, 'deg'), [-1 0 0; 0 -1 0; 0 0 1],'absTol',1e-10);
    
    syms q
    R = rotz(q);
    verifyInstanceOf(tc, R, 'sym');
    verifySize(tc,R, [3 3]);
    tc.verifyEqual(simplify(det(R)), sym(1));
    
    %test for non-scalar input
    verifyError(tc, @()rotz([1 2 3]),'SMTB:rotz:badarg');
end

function trotx_test(tc)
    tc.verifyEqual(trotx(0), eye(4,4),'absTol',1e-10);
    tc.verifyEqual(trotx(pi/2), [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1],'absTol',1e-10);
    tc.verifyEqual(trotx(pi), [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1],'absTol',1e-10);
    
    tc.verifyEqual(trotx(90, 'deg'), [1 0 0 0; 0 0 -1 0; 0 1 0 0; 0 0 0 1],'absTol',1e-10);
    tc.verifyEqual(trotx(180, 'deg'), [1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 1],'absTol',1e-10);
    
    %test for non-scalar input
    verifyError(tc, @()trotx([1 2 3; 0 0 0 1]),'MATLAB:catenate:dimensionMismatch');
end

function troty_test(tc)
    tc.verifyEqual(troty(0), eye(4,4),'absTol',1e-10);
    tc.verifyEqual(troty(pi/2), [0 0 1 0; 0 1 0 0; -1 0 0 0; 0 0 0 1],'absTol',1e-10);
    tc.verifyEqual(troty(pi), [-1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1],'absTol',1e-10);
    
    tc.verifyEqual(troty(90, 'deg'), [0 0 1 0; 0 1 0 0; -1 0 0 0; 0 0 0 1],'absTol',1e-10);
    tc.verifyEqual(troty(180, 'deg'), [-1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1],'absTol',1e-10);
    %test for non-scalar input
    verifyError(tc, @()troty([1 2 3; 0 0 0 1]),'MATLAB:catenate:dimensionMismatch');
end

function trotz_test(tc)
    tc.verifyEqual(trotz(0), eye(4,4),'absTol',1e-10);
    tc.verifyEqual(trotz(pi/2), [0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1],'absTol',1e-10);
    tc.verifyEqual(trotz(pi), [-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1],'absTol',1e-10);
    
    tc.verifyEqual(trotz(90, 'deg'), [0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1],'absTol',1e-10);
    tc.verifyEqual(trotz(180, 'deg'), [-1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1],'absTol',1e-10);
    %test for non-scalar input
    verifyError(tc, @()trotz([1 2 3; 0 0 0 1]),'MATLAB:catenate:dimensionMismatch');
end

function transl_test(tc)
    
    % transl(P) -> T
    tc.verifyEqual(transl(1, 2, 3), [1 0 0 1; 0 1 0 2; 0 0 1 3; 0 0 0 1], 'AbsTol', 1e-10);
    tc.verifyEqual(transl([1, 2, 3]), [1 0 0 1; 0 1 0 2; 0 0 1 3; 0 0 0 1], 'AbsTol', 1e-10);
    
    x = rand(4,3);
    T = transl(x);
    tc.verifyEqual(T(1:3,4,1), x(1,:)', 'AbsTol', 1e-10);
    tc.verifyEqual(T(1:3,4,4), x(4,:)', 'AbsTol', 1e-10);
    
    % transl(T) -> P
    tc.verifyEqual(transl([1 0 0 1; 0 1 0 2; 0 0 1 3; 0 0 0 1]), [1 2 3]', 'AbsTol', 1e-10);
    [a,b,c] = transl([1 0 0 1; 0 1 0 2; 0 0 1 3; 0 0 0 1]);
    tc.verifyEqual(a, 1);
    tc.verifyEqual(b, 2);
    tc.verifyEqual(c, 3);
    
    tc.verifyEqual(transl(T), x, 'AbsTol', 1e-10);
end




%% angle/vector form
%    angvec2r                   - angle/vector to RM
function angvec2r_test(tc)
    
    tc.verifyEqual(angvec2r( pi/2, [1 0 0]), rotx(pi/2),'absTol',1e-10);
    tc.verifyEqual(angvec2r( pi/2, [0 1 0]), roty(pi/2),'absTol',1e-10);
    tc.verifyEqual(angvec2r( pi/2, [0 0 1]), rotz(pi/2),'absTol',1e-10);
    
    tc.verifyEqual(angvec2r(0, [1 0 0]), eye(3,3),'absTol',1e-10);
    tc.verifyEqual(angvec2r(0, [0 1 0]), eye(3,3),'absTol',1e-10);
    tc.verifyEqual(angvec2r(0, [0 0 1]), eye(3,3),'absTol',1e-10);
    tc.verifyEqual(angvec2r(0, [0 0 0]), eye(3,3),'absTol',1e-10);
    
    verifyError(tc, @()angvec2r(1, [0 0 0]),'SMTB:angvec2r:badarg');
    
    verifyError(tc, @()angvec2r([1,2,3],0.1),'SMTB:angvec2r:badarg');
    verifyError(tc, @()angvec2r(1),'SMTB:angvec2r:badarg');
end

%    angvec2tr                  - angle/vector to HT
function angvec2tr_test(tc)
    
    tc.verifyEqual(angvec2tr( pi/2, [1 0 0]), trotx(pi/2),'absTol',1e-10);
    tc.verifyEqual(angvec2tr( pi/2, [0 1 0]), troty(pi/2),'absTol',1e-10);
    tc.verifyEqual(angvec2tr( pi/2, [0 0 1]), trotz(pi/2),'absTol',1e-10);
    
    tc.verifyEqual(angvec2tr(0, [1 0 0]), eye(4,4),'absTol',1e-10);
    tc.verifyEqual(angvec2tr(0, [0 1 0]), eye(4,4),'absTol',1e-10);
    tc.verifyEqual(angvec2tr(0, [0 0 1]), eye(4,4),'absTol',1e-10);
    tc.verifyEqual(angvec2tr(0, [0 0 0]), eye(4,4),'absTol',1e-10);
    
    verifyError(tc, @()angvec2tr(1, [0 0 0]),'SMTB:angvec2r:badarg');
    
    verifyError(tc, @()angvec2tr([1,2,3],0.1),'SMTB:angvec2r:badarg');
    verifyError(tc, @()angvec2tr(1),'SMTB:angvec2tr:badarg');
end

%    tr2angvec                  - HT/RM to angle/vector form
function tr2angvec_test(tc)
    % null rotation
    % - vector isn't defined here, but RTB sets it (0 0 0)
    [theta, v] = tr2angvec(eye(3,3));
    tc.verifyEqual(theta, 0.0, 'absTol',1e-6);
    tc.verifyEqual(v, [0 0 0], 'absTol',1e-6);
    
    tr2angvec(eye(3,3))
    
    % canonic rotations
    [theta, v] = tr2angvec(rotx(pi/2));
    tc.verifyEqual(theta, pi/2, 'absTol',1e-6);
    tc.verifyEqual(v, [1 0 0], 'absTol',1e-6);
    
    [theta, v] = tr2angvec(roty(pi/2));
    tc.verifyEqual(theta, pi/2, 'absTol',1e-6);
    tc.verifyEqual(v, [0 1 0], 'absTol',1e-6);
    
    [theta, v] = tr2angvec(rotz(pi/2));
    tc.verifyEqual(theta, pi/2, 'absTol',1e-6);
    tc.verifyEqual(v, [0 0 1], 'absTol',1e-6);
    
    % null rotation
    [theta, v] = tr2angvec(eye(4,4));
    tc.verifyEqual(theta, 0.0, 'absTol',1e-6);
    tc.verifyEqual(v, [0 0 0], 'absTol',1e-6);
    
    % canonic rotations
    [theta, v] = tr2angvec(trotx(pi/2));
    tc.verifyEqual(theta, pi/2, 'absTol',1e-6);
    tc.verifyEqual(v, [1 0 0], 'absTol',1e-6);
    
    [theta, v] = tr2angvec(troty(pi/2));
    tc.verifyEqual(theta, pi/2, 'absTol',1e-6);
    tc.verifyEqual(v, [0 1 0], 'absTol',1e-6);
    
    [theta, v] = tr2angvec(trotz(pi/2));
    tc.verifyEqual(theta, pi/2, 'absTol',1e-6);
    tc.verifyEqual(v, [0 0 1], 'absTol',1e-6);
    
    [theta, v] = tr2angvec(roty(pi/2), 'deg');
    tc.verifyEqual(theta, 90, 'absTol',1e-6);
    tc.verifyEqual(v, [0 1 0], 'absTol',1e-6);
    
    R = cat(3, rotx(pi/2), roty(pi/2), rotz(pi/2));
    [theta, v] = tr2angvec(R);
    tc.verifyEqual(theta, pi/2*[1 1 1]', 'absTol',1e-6);
    tc.verifyEqual(v, eye(3,3), 'absTol',1e-6);
    
    T = cat(3, trotx(pi/2), troty(pi/2), trotz(pi/2));
    [theta, v] = tr2angvec(T);
    tc.verifyEqual(theta, pi/2*[1 1 1]', 'absTol',1e-6);
    tc.verifyEqual(v, eye(3,3), 'absTol',1e-6);
    
    %test for scalar input
    verifyError(tc, @()tr2angvec(1), 'SMTB:tr2angvec:badarg');
end



%% 3-angle forms

function eul2r_test(tc)
    
    % ZYZ
    r2d = 180/pi;
    
    R = rotz(0.1) * roty(0.2) * rotz(0.3);
    
    tc.verifyEqual(eul2r(0.1, 0.2, 0.3), R, 'absTol',1e-10);
    tc.verifyEqual(eul2r([0.1, 0.2, 0.3]), R, 'absTol',1e-10);
    tc.verifyEqual(eul2r(0.1*r2d, 0.2*r2d, 0.3*r2d, 'deg'), R, 'absTol',1e-10);
    tc.verifyEqual(eul2r([0.1, 0.2, 0.3]*r2d, 'deg'), R, 'absTol',1e-10);
    
    % trajectory case
    Rs = eul2r( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]);
    verifySize(tc, Rs, [3 3 3]);
    tc.verifyEqual(Rs(:,:,2), R, 'absTol',1e-10);
    
    Rs = eul2r( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]*r2d, 'deg');
    verifySize(tc, Rs, [3 3 3]);
    tc.verifyEqual(Rs(:,:,2), R, 'absTol',1e-10);
    
    %test for scalar input
    verifyError(tc, @()eul2r(1),'SMTB:eul2r:badarg');
end

%    eul2tr                     - Euler angles to HT
function eul2tr_test(tc)
    r2d = 180/pi;
    
    T = trotz(0.1) * troty(0.2) * trotz(0.3);
    
    tc.verifyEqual(eul2tr(0.1, 0.2, 0.3), T, 'absTol',1e-10);
    tc.verifyEqual(eul2tr([0.1, 0.2, 0.3]), T, 'absTol',1e-10);
    tc.verifyEqual(eul2tr(0.1*r2d, 0.2*r2d, 0.3*r2d, 'deg'), T, 'absTol',1e-10);
    tc.verifyEqual(eul2tr([0.1, 0.2, 0.3]*r2d, 'deg'), T, 'absTol',1e-10);
    
    % trajectory case
    Ts = eul2tr( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]);
    verifySize(tc, Ts, [4 4 3]);
    tc.verifyEqual(Ts(:,:,2), T, 'absTol',1e-10);
    
    Ts = eul2tr( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]*r2d, 'deg');
    verifySize(tc, Ts, [4 4 3]);
    tc.verifyEqual(Ts(:,:,2), T, 'absTol',1e-10);
    
    %test for scalar input
    verifyError(tc, @()eul2tr(1),'SMTB:eul2r:badarg');
end

function rpy2r_test(tc)
    
    r2d = 180/pi;
    
    %% default zyx order
    R = rotz(0.3) * roty(0.2) * rotx(0.1);
    
    tc.verifyEqual(rpy2r(0.1, 0.2, 0.3), R, 'absTol',1e-10);
    tc.verifyEqual(rpy2r([0.1, 0.2, 0.3]), R, 'absTol',1e-10);
    tc.verifyEqual(rpy2r(0.1*r2d, 0.2*r2d, 0.3*r2d, 'deg'), R, 'absTol',1e-10);
    tc.verifyEqual(rpy2r([0.1, 0.2, 0.3]*r2d, 'deg'), R, 'absTol',1e-10);
    
    % trajectory case
    Rs = rpy2r( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]);
    verifySize(tc, Rs, [3 3 3]);
    tc.verifyEqual(Rs(:,:,2), R, 'absTol',1e-10);
    
    Rs = rpy2r( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]*r2d, 'deg');
    verifySize(tc, Rs, [3 3 3]);
    tc.verifyEqual(Rs(:,:,2), R, 'absTol',1e-10);
    
    %% xyz order
    
    R = rotx(0.3) * roty(0.2) * rotz(0.1);
    
    tc.verifyEqual(rpy2r(0.1, 0.2, 0.3, 'xyz'), R, 'absTol',1e-10);
    tc.verifyEqual(rpy2r([0.1, 0.2, 0.3], 'xyz'), R, 'absTol',1e-10);
    tc.verifyEqual(rpy2r(0.1*r2d, 0.2*r2d, 0.3*r2d, 'deg', 'xyz'), R, 'absTol',1e-10);
    tc.verifyEqual(rpy2r([0.1, 0.2, 0.3]*r2d, 'deg', 'xyz'), R, 'absTol',1e-10);
    
    % trajectory case
    Rs = rpy2r( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3], 'xyz');
    verifySize(tc, Rs, [3 3 3]);
    tc.verifyEqual(Rs(:,:,2), R, 'absTol',1e-10);
    
    Rs = rpy2r( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]*r2d, 'xyz', 'deg');
    verifySize(tc, Rs, [3 3 3]);
    tc.verifyEqual(Rs(:,:,2), R, 'absTol',1e-10);
    
    %% yxz order
    
    R = roty(0.3) * rotx(0.2) * rotz(0.1);
    
    tc.verifyEqual(rpy2r(0.1, 0.2, 0.3, 'yxz'), R, 'absTol',1e-10);
    tc.verifyEqual(rpy2r([0.1, 0.2, 0.3], 'yxz'), R, 'absTol',1e-10);
    tc.verifyEqual(rpy2r(0.1*r2d, 0.2*r2d, 0.3*r2d, 'deg', 'yxz'), R, 'absTol',1e-10);
    tc.verifyEqual(rpy2r([0.1, 0.2, 0.3]*r2d, 'deg', 'yxz'), R, 'absTol',1e-10);
    
    % trajectory case
    Rs = rpy2r( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3], 'yxz');
    verifySize(tc, Rs, [3 3 3]);
    tc.verifyEqual(Rs(:,:,2), R, 'absTol',1e-10);
    
    Rs = rpy2r( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]*r2d, 'yxz', 'deg');
    verifySize(tc, Rs, [3 3 3]);
    tc.verifyEqual(Rs(:,:,2), R, 'absTol',1e-10);
    
    %test for scalar input
    verifyError(tc, @()rpy2tr(1),'SMTB:rpy2r:badarg');
end

function rpy2tr_test(tc)
    
    r2d = 180/pi;
    
    T = trotz(0.3) * troty(0.2) * trotx(0.1);
    seq = 'zyx';
    tc.verifyEqual(rpy2tr(0.1, 0.2, 0.3, seq), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr([0.1, 0.2, 0.3], seq), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr(0.1*r2d, 0.2*r2d, 0.3*r2d, 'deg', seq), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr([0.1, 0.2, 0.3]*r2d, 'deg', seq), T, 'absTol',1e-10);
    
    
    T = trotx(0.3) * troty(0.2) * trotz(0.1);
    seq = 'xyz';
    tc.verifyEqual(rpy2tr(0.1, 0.2, 0.3, seq), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr([0.1, 0.2, 0.3], seq), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr(0.1*r2d, 0.2*r2d, 0.3*r2d, 'deg', seq), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr([0.1, 0.2, 0.3]*r2d, 'deg', seq), T, 'absTol',1e-10);
    
    T = troty(0.3) * trotx(0.2) * trotz(0.1);
    seq = 'yxz';
    tc.verifyEqual(rpy2tr(0.1, 0.2, 0.3, seq), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr([0.1, 0.2, 0.3], seq), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr(0.1*r2d, 0.2*r2d, 0.3*r2d, 'deg', seq), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr([0.1, 0.2, 0.3]*r2d, 'deg', seq), T, 'absTol',1e-10);
    
    
    % trajectory case
    T = trotz(0.3) * troty(0.2) * trotx(0.1);
    
    Ts = rpy2tr( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]);
    verifySize(tc, Ts, [4 4 3]);
    tc.verifyEqual(Ts(:,:,2), T, 'absTol',1e-10);
    
    Ts = rpy2tr( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]*r2d, 'deg');
    verifySize(tc, Ts, [4 4 3]);
    tc.verifyEqual(Ts(:,:,2), T, 'absTol',1e-10);
    
    T = trotx(0.3) * troty(0.2) * trotz(0.1);
    
    tc.verifyEqual(rpy2tr(0.1, 0.2, 0.3, 'xyz'), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr([0.1, 0.2, 0.3], 'xyz'), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr(0.1*r2d, 0.2*r2d, 0.3*r2d, 'deg', 'xyz'), T, 'absTol',1e-10);
    tc.verifyEqual(rpy2tr([0.1, 0.2, 0.3]*r2d, 'deg', 'xyz'), T, 'absTol',1e-10);
    
    % trajectory case
    Ts = rpy2tr( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3], 'xyz');
    verifySize(tc, Ts, [4 4 3]);
    tc.verifyEqual(Ts(:,:,2), T, 'absTol',1e-10);
    
    Ts = rpy2tr( [0.1, 0.2, 0.3; 0.1 0.2 0.3; 0.1 0.2 0.3]*r2d, 'xyz', 'deg');
    verifySize(tc, Ts, [4 4 3]);
    tc.verifyEqual(Ts(:,:,2), T, 'absTol',1e-10);
    
    %test for scalar input
    verifyError(tc, @()rpy2tr(1),'SMTB:rpy2r:badarg');
end

function tr2eul_test(tc)
    
    eul = [0.1 0.2 0.3];
    R = eul2r(eul);
    tc.verifyEqual(tr2eul(R), eul,'absTol',1e-10);
    tc.verifyEqual(tr2eul(R, 'deg'), eul*180/pi,'absTol',1e-10);
    
    Rs = cat(3, R, R, R, R);
    x = tr2eul(Rs);
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), eul,'absTol',1e-10);
    x = tr2eul(Rs, 'deg');
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), eul*180/pi,'absTol',1e-10);
    
    T = eul2tr(eul);
    tc.verifyEqual(tr2eul(T), eul,'absTol',1e-10);
    tc.verifyEqual(tr2eul(T, 'deg'), eul*180/pi,'absTol',1e-10);
    
    Ts = cat(3, T, T, T, T);
    x = tr2eul(Ts);
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), eul,'absTol',1e-10);
    x = tr2eul(Ts, 'deg');
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), eul*180/pi,'absTol',1e-10);
    
    % test singularity case
    eul = [0.1 0 0.3];
    R = eul2r(eul);
    tc.verifyEqual(eul2r( tr2eul(R) ), R,'absTol',1e-10);
    tc.verifyEqual(eul2r( tr2eul(R, 'deg'), 'deg'), R,'absTol',1e-10);
    
    %test for scalar input
    verifyError(tc, @()tr2eul(1),'SMTB:tr2eul:badarg');

    % test flip
    eul = [-0.1 0.2 0.3];
    R = eul2r(eul);
    eul2 = tr2eul(R, 'flip');
    tc.verifyTrue(eul2(1) > 0);
    tc.verifyEqual(eul2r(eul2), R,'absTol',1e-10);
end

function tr2rpy_test(tc)
    rpy = [0.1 0.2 0.3];
    R = rpy2r(rpy);
    tc.verifyEqual(tr2rpy(R), rpy,'absTol',1e-10);
    tc.verifyEqual(tr2rpy(R, 'deg'), rpy*180/pi,'absTol',1e-10);
    
    Rs = cat(3, R, R, R, R);
    x = tr2rpy(Rs);
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), rpy,'absTol',1e-10);
    x = tr2rpy(Rs, 'deg');
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), rpy*180/pi,'absTol',1e-10);
    
    T = rpy2tr(rpy);
    tc.verifyEqual(tr2rpy(T), rpy,'absTol',1e-10);
    tc.verifyEqual(tr2rpy(T, 'deg'), rpy*180/pi,'absTol',1e-10);
    
    Ts = cat(3, T, T, T, T);
    x = tr2rpy(Ts);
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), rpy,'absTol',1e-10);
    x = tr2rpy(Ts, 'deg');
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), rpy*180/pi,'absTol',1e-10);
    
    % xyz order
    R = rpy2r(rpy, 'xyz');
    tc.verifyEqual(tr2rpy(R, 'xyz'), rpy,'absTol',1e-10);
    tc.verifyEqual(tr2rpy(R, 'deg', 'xyz'), rpy*180/pi,'absTol',1e-10);
    
    Rs = cat(3, R, R, R, R);
    x = tr2rpy(Rs, 'xyz');
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), rpy,'absTol',1e-10);
    x = tr2rpy(Rs, 'deg', 'xyz');
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), rpy*180/pi,'absTol',1e-10);
    
    T = rpy2tr(rpy, 'xyz');
    tc.verifyEqual(tr2rpy(T, 'xyz'), rpy,'absTol',1e-10);
    tc.verifyEqual(tr2rpy(T, 'deg', 'xyz'), rpy*180/pi,'absTol',1e-10);
    
    Ts = cat(3, T, T, T, T);
    x = tr2rpy(Ts, 'xyz');
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), rpy,'absTol',1e-10);
    x = tr2rpy(Ts, 'deg', 'xyz');
    verifySize(tc, x, [4 3]);
    tc.verifyEqual(x(2,:), rpy*180/pi,'absTol',1e-10);
    
    % corner cases
    seq = 'zyx';
    ang = [pi 0 0];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 pi 0];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 0 pi];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 pi/2 0]; % singularity
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 -pi/2 0];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    
    seq = 'xyz';
    ang = [pi 0 0];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 pi 0];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 0 pi];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 pi/2 0]; % singularity
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 -pi/2 0];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    
    seq = 'yxz';
    ang = [pi 0 0];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 pi 0];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 0 pi];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 pi/2 0]; % singularity
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    ang = [0 -pi/2 0];
    a = rpy2tr(ang, seq);
    tc.verifyEqual(rpy2tr(tr2rpy(a, seq), seq), a, 'absTol',1e-10);
    
    %test for scalar input
    verifyError(tc, @()tr2rpy(1),'SMTB:tr2rpy:badarg');
end



%    oa2r                       - orientation and approach vector to RM
function oa2r_test(tc)
    %Unit test for oa2r with variables ([0 1 0] & [0 0 1])
    tc.verifyEqual(oa2r([0 1 0], [0 0 1]),...
        [1     0     0
        0     1     0
        0     0     1],'absTol',1e-10);
    %test for scalar input
    verifyError(tc, @()oa2r(1),'SMTB:oa2r:badarg');
end

%    oa2tr                      - orientation and approach vector to HT
function oa2tr_test(tc)
    %Unit test for oa2tr with variables ([0 1 0] & [0 0 1])
    tc.verifyEqual(oa2tr([0 1 0], [0 0 1]),...
        [1     0     0     0
        0     1     0     0
        0     0     1     0
        0     0     0     1],'absTol',1e-10);
    %test for scalar input
    verifyError(tc, @()oa2tr(1),'SMTB:oa2tr:badarg');
end


function trchain_test(tc)
    a1 = 0;
    
    T = trchain('Tx(a1) Ty(a1) Tz(a1) Rx(a1) Ry(a1) Rz(a1)');
    tc.verifyEqual(T, eye(4,4), 'abstol', 1e-10);
    
    a1 = 1; a2 = 2; a3 = 3;
    tc.verifyEqual( trchain('Tx(a1) Ty(a2) Tz(a3)'), transl(1,2,3), 'abstol', 1e-10);
    
    a1 = 0.3; a2 = 0.4; a3 = 0.5;
    tc.verifyEqual( trchain('Rx(a1) Ry(a2) Rz(a3)'), trotx(0.3)*troty(0.4)*trotz(0.5), 'abstol', 1e-10);
    
    tc.verifyEqual( trchain('Rx(q1) Ry(q2) Rz(q3)', [.3,.4,.5]), trotx(0.3)*troty(0.4)*trotz(0.5), 'abstol', 1e-10);
    
    syms q1 q2 q3 a1 a2 a3
    tc.verifyEqual( trchain('Rx(q1) Tx(a1) Ry(q2) Ty(a2) Rz(q3) Tz(a3)', [q1 q2 q3]), trotx(q1)*transl(a1,0,0)*troty(q2)*transl(0,a2,0)*trotz(q3)*transl(0,0,a3) );

    syms q1(t) q2(t) q3(t) t a1 a2 a3
    tc.verifyEqual( trchain('Rx(q1) Tx(a1) Ry(q2) Ty(a2) Rz(q3) Tz(a3)', [q1 q2 q3]), formula(trotx(q1)*transl(a1,0,0)*troty(q2)*transl(0,a2,0)*trotz(q3)*transl(0,0,a3)) );
end


function trchain2_test(tc)
    a1 = 0;
    
    T = trchain2('Tx(a1) Ty(a1) R(a1)');
    tc.verifyEqual(T, eye(3,3), 'abstol', 1e-10);
    
    a1 = 1; a2 = 2;
    tc.verifyEqual( trchain2('Tx(a1) Ty(a2)'), transl2(1,2), 'abstol', 1e-10);
    
    a1 = 0.3; a2 = 0.4;
    % R() is the same as Rz()
    tc.verifyEqual( trchain2('R(a1) Rz(a2)'), trot2(0.3)*trot2(0.4), 'abstol', 1e-10);
    
    syms q1 a1 a2
    tc.verifyEqual( trchain2('R(q1) Tx(a1) Ty(a2)', [q1]), trot2(q1)*transl2(a1,0)*transl2(0,a2) );
end

function trinterp_test(tc)
    %% between two transforms
    T0 = transl(1,2,3);
    T1 = transl(-1,-2,-3)*trotx(pi);
    Tm = trotx(pi/2);
    
    T = trinterp(T0, T1, 0);
    tc.verifyEqual(size(T), [4 4]);
    tc.verifyEqual(T, T0, 'abstol', 1e-10);
    
    tc.verifyEqual(trinterp(T0, T1, 1), T1, 'abstol', 1e-10);
    tc.verifyEqual(trinterp(T0, T1, 0.5), Tm, 'abstol', 1e-10);
    
    T = trinterp(T0, T1, [0.5 0 1]);
    tc.verifyEqual(size(T), [4 4 3]);
    tc.verifyEqual(T(:,:,1), Tm, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,2), T0, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,3), T1, 'abstol', 1e-10);
    
    T = trinterp(T0, T1, 3);   % interpolate in 3 steps
    tc.verifyEqual(size(T), [4 4 3]);
    tc.verifyEqual(T(:,:,1), T0, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,2), Tm, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,3), T1, 'abstol', 1e-10);
    
    %% between identity and transform
    T0 = eye(4,4);
    T1 = transl(2,4,6)*trotx(pi);
    Tm = transl(1,2,3)*trotx(pi/2);
    
    T = trinterp(T1, 0);
    tc.verifyEqual(size(T), [4 4]);
    tc.verifyEqual(T, T0, 'abstol', 1e-10);
    
    tc.verifyEqual(trinterp(T1, 1), T1, 'abstol', 1e-10);
    tc.verifyEqual(trinterp(T1, 0.5), Tm, 'abstol', 1e-10);
    
    T = trinterp(T1, [0.5 0 1]);
    tc.verifyEqual(size(T), [4 4 3]);
    tc.verifyEqual(T(:,:,1), Tm, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,2), T0, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,3), T1, 'abstol', 1e-10);
    
    T = trinterp(T1, 3);   % interpolate in 3 steps
    tc.verifyEqual(size(T), [4 4 3]);
    tc.verifyEqual(T(:,:,1), T0, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,2), Tm, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,3), T1, 'abstol', 1e-10);
    
    tc.verifyError( @() trinterp(T0, T1, -1), 'SMTB:trinterp:badarg');
    tc.verifyError( @() trinterp(T0, T1, 1.7), 'SMTB:trinterp:badarg');
    tc.verifyError( @() trinterp(T0), 'SMTB:trinterp:badarg');

end


function trinterp2_test(tc)
    %% between two transforms
    T0 = transl2(1,2);
    T1 = transl2(-1,-2)*trot2(pi);
    Tm = trot2(pi/2);
    
    T = trinterp2(T0, T1, 0);
    tc.verifyEqual(size(T), [3 3]);
    tc.verifyEqual(T, T0, 'abstol', 1e-10);
    
    tc.verifyEqual(trinterp2(T0, T1, 1), T1, 'abstol', 1e-10);
    tc.verifyEqual(trinterp2(T0, T1, 0.5), Tm, 'abstol', 1e-10);
    
    T = trinterp2(T0, T1, [0.5 0 1]);
    tc.verifyEqual(size(T), [3 3 3]);
    tc.verifyEqual(T(:,:,1), Tm, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,2), T0, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,3), T1, 'abstol', 1e-10);
    
    T = trinterp2(T0, T1, 3);
    tc.verifyEqual(size(T), [3 3 3]);
    tc.verifyEqual(T(:,:,1), T0, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,2), Tm, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,3), T1, 'abstol', 1e-10);
    
    %% between identity and transform
    T0 = eye(3, 3);
    T1 = transl2(2,4)*trot2(pi);
    Tm = transl2(1,2)*trot2(pi/2);
    
    T = trinterp2(T1, 0);
    tc.verifyEqual(size(T), [3 3]);
    tc.verifyEqual(T, T0, 'abstol', 1e-10);
    
    tc.verifyEqual(trinterp2(T0, T1, 1), T1, 'abstol', 1e-10);
    tc.verifyEqual(trinterp2(T0, T1, 0.5), Tm, 'abstol', 1e-10);
    
    T = trinterp2(T0, T1, [0.5 0 1]);
    tc.verifyEqual(size(T), [3 3 3]);
    tc.verifyEqual(T(:,:,1), Tm, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,2), T0, 'abstol', 1e-10);
    tc.verifyEqual(T(:,:,3), T1, 'abstol', 1e-10);
end


%    trnorm                     - normalize HT
function trnorm_test(tc)
    
    R = [0.9 0 0; .2 .6 .3; .1 .2 .4]';
    tc.verifyEqual(det(trnorm(R)), 1, 'absTol', 1e-14);
    
    t = [1 2 3]';
    T = rt2tr(R, t);
    Tn = trnorm(T);
    tc.verifyEqual(det(trnorm(t2r(Tn))), 1, 'absTol', 1e-14);
    tc.verifyEqual(Tn(1:3,4), t);
    
    % vector input
    RR = cat(3, R, R, R, R);
    RRn = trnorm(RR);
    verifySize(tc, RRn, [3 3 4]);
    tc.verifyEqual(det(RRn(:,:,1)), 1, 'absTol', 1e-14);
    tc.verifyEqual(det(RRn(:,:,1)), 1, 'absTol', 1e-14);
    
    %HACK)    tc.verifyEqual(arrayfun( @(x) det(trnorm(x)), RR
    
    %test for scalar input
    verifyError(tc, @()trnorm(1),'SMTB:trnorm:badarg');
end

function trprint_test(tc)
    
    % null case
    
    trprint(eye(4,4))
    
    a = transl([1,2,3]) * eul2tr([.1, .2, .3]);
    
    trprint(a);

    trprint('a') % equivalent to trprint a on the command line
    
    s = evalc( 'trprint(a)' );
    tc.verifyClass(s, 'char');
    tc.verifyEqual(size(s,1), 1);
    s = evalc( 'trprint(cat(3, a, a, a))' );
    tc.verifyClass(s, 'char');
    tc.verifyEqual(size(s,1), 1);
    
    s = trprint(a);
    tc.verifyClass(s, 'char');
    tc.verifyEqual(size(s,1), 1);
    
    tc.verifyClass(s, 'char');
    tc.verifyEqual(size(s,1), 1);
    
    trprint(a, 'euler');
    trprint(a, 'euler', 'radian');
    trprint(a, 'rpy');
    trprint(a, 'rpy', 'radian');
    trprint(a, 'rpy', 'radian', 'xyz');
    trprint(a, 'rpy', 'radian', 'zyx');
    trprint(a, 'rpy', 'radian', 'yxz');
    
    trprint(a, 'angvec');
    trprint(a, 'angvec', 'radian');
    trprint(a, 'angvec', 'radian', 'fmt', '%g');
    trprint(a, 'angvec', 'radian', 'fmt', '%g', 'label', 'bob');
    
    % vector case
    
    a = cat(3, a, a, a);
    trprint(a);
    
    s = evalc( 'trprint(a)' );
    tc.verifyTrue(isa(s, 'char') );
    tc.verifyEqual( length(regexp(s, '\n', 'match')), 4);
    
    trprint(a, 'euler');
    trprint(a, 'euler', 'radian');
    trprint(a, 'rpy');
    trprint(a, 'rpy', 'radian');
    trprint(a, 'rpy', 'radian', 'xyz');
    trprint(a, 'rpy', 'radian', 'zyx');
    trprint(a, 'rpy', 'radian', 'yxz');
    
    trprint(a, 'angvec');
    trprint(a, 'angvec', 'radian');
    trprint(a, 'angvec', 'radian', 'fmt', '%g');
    trprint(a, 'angvec', 'radian', 'fmt', '%g', 'label', 'bob');
    
    % write to a file
    f = fopen('test.txt', 'w');
    trprint(a, 'fid', f);
    fclose(f);
    % read it back
    f = fopen('test.txt', 'r');
    s2 = fread(f, '*char');
    fclose(f);
    delete 'test.txt';
    tc.verifyTrue( all(s(:) == s2) );
    
end


function trscale_test(tc)
    
    tc.verifyEqual( trscale(1), eye(4,4) );
    tc.verifyEqual( trscale([1 2 3]), diag([1 2 3 1]) );
    tc.verifyEqual( trscale(1, 2, 3), diag([1 2 3 1]) );
    
end

function vex_test(tc)
    S = [
        0    -3     2
        3     0    -1
        -2     1     0
        ];
    
    tc.verifyEqual( vex(S), [1 2 3]');
    
    tc.verifyEqual( vex(-S), -[1  2 3]');
end


function skew_test(tc)
    R = skew([1 2 3]);
    
    tc.verifyTrue( isrot(R) );  % check size
    
    tc.verifyEqual( norm(R'+ R), 0, 'absTol', 1e-10); % check is skew
    
    tc.verifyEqual( vex(R), [1 2 3]'); % check contents, vex already verified

    tc.verifyError( @() skew([1 2]), 'SMTB:skew:badarg')
end

function vexa_test(tc)
    S = [
        0    -6     5     1
        6     0    -4     2
        -5     4     0     3
        0     0     0     0
        ];
    
    tc.verifyEqual( vexa(S), [1:6]');
    
    S = [
        0     6     5     1
        -6     0     4    -2
        -5    -4     0     3
        0     0     0     0
        ];
    tc.verifyEqual( vexa(S), [1 -2 3 -4 5 -6]');
end


function skewa_test(tc)
    T = skewa([3 4 5]);
    
    tc.verifyTrue( ishomog2(T) );  % check size
    
    R = T(1:2,1:2);
    tc.verifyEqual( norm(R'+ R), 0, 'absTol', 1e-10); % check is skew
    
    tc.verifyEqual( vexa(T), [3 4 5]'); % check contents, vexa already verified
    tc.verifyError( @() skewa([1 2]), 'SMTB:skewa:badarg')

end

function trlog_test(tc)
    %unit tests for matrix expon stuff
    
    %%% SO(3) tests
    % zero rotation case
    tc.verifyEqual(trlog( eye(3,3) ), skew([0 0 0]), 'absTol', 1e-6);
    
    % rotation by pi case
    tc.verifyEqual(trlog( rotx(pi) ), skew([pi 0 0]), 'absTol', 1e-6);
    tc.verifyEqual(trlog( roty(pi) ), skew([0 pi 0]), 'absTol', 1e-6);
    tc.verifyEqual(trlog( rotz(pi) ), skew([0 0 pi]), 'absTol', 1e-6);
    
    % general case
    tc.verifyEqual(trlog( rotx(0.2) ), skew([0.2 0 0]), 'absTol', 1e-6);
    tc.verifyEqual(trlog( roty(0.3) ), skew([0 0.3 0]), 'absTol', 1e-6);
    tc.verifyEqual(trlog( rotz(0.4) ), skew([0 0 0.4]), 'absTol', 1e-6);
    
    
    R = rotx(0.2) * roty(0.3) * rotz(0.4);
    [th,w] = trlog(R);
    tc.verifyEqual( logm(R), skew(th*w), 'absTol', 1e-10)
    
    %%% SE(3) tests
    
    % pure translation
    tc.verifyEqual(trlog( transl([1 2 3]) ), ...
        [0 0 0 1; 0 0 0 2; 0 0 0 3; 0 0 0 0], 'absTol', 1e-6);
    
    % mixture
    T = transl([1 2 3])*trotx(0.3);
    tc.verifyEqual(trlog(T), logm(T), 'absTol', 1e-6);
    
    T = transl([1 2 3])*troty(0.3);
    tc.verifyEqual(trlog(T), logm(T), 'absTol', 1e-6);
    
    [th,w] = trlog(T);
    tc.verifyEqual( logm(T), skewa(th*w), 'absTol', 1e-10)
    
    
    verifyError(tc, @()trlog(0),'SMTB:trlog:badarg');
end


function trexp_test(tc)
    %unit tests for matrix log stuff
    
    %%% SO(3) tests
    
    %% so(3)
    
    % zero rotation case
    tc.verifyEqual(trexp(skew([0 0 0])), eye(3,3), 'absTol', 1e-6);
    
    %% so(3), theta
    
    tc.verifyEqual(trexp(skew([0 0 0]), 1), eye(3,3), 'absTol', 1e-6);
    
    % rotation by pi case
    tc.verifyEqual(trexp(skew([pi 0 0])), rotx(pi), 'absTol', 1e-6);
    tc.verifyEqual(trexp(skew([0 pi 0])), roty(pi), 'absTol', 1e-6);
    tc.verifyEqual(trexp(skew([0 0 pi])), rotz(pi), 'absTol', 1e-6);
    
    % general case
    tc.verifyEqual(trexp(skew([0.2 0 0])), rotx(0.2), 'absTol', 1e-6);
    tc.verifyEqual(trexp(skew([0 0.3 0])), roty(0.3), 'absTol', 1e-6);
    tc.verifyEqual(trexp(skew([0 0 0.4])), rotz(0.4), 'absTol', 1e-6);
    
    tc.verifyEqual(trexp(skew([1 0 0]), 0.2), rotx(0.2), 'absTol', 1e-6);
    tc.verifyEqual(trexp(skew([0 1 0]), 0.3), roty(0.3), 'absTol', 1e-6);
    tc.verifyEqual(trexp(skew([0 0 1]), 0.4), rotz(0.4), 'absTol', 1e-6);
    
    tc.verifyEqual(trexp([1 0 0], 0.2), rotx(0.2), 'absTol', 1e-6);
    tc.verifyEqual(trexp([0 1 0], 0.3), roty(0.3), 'absTol', 1e-6);
    tc.verifyEqual(trexp([0 0 1], 0.4), rotz(0.4), 'absTol', 1e-6);
    
    tc.verifyEqual(trexp([1 0 0]*0.2), rotx(0.2), 'absTol', 1e-6);
    tc.verifyEqual(trexp([0 1 0]*0.3), roty(0.3), 'absTol', 1e-6);
    tc.verifyEqual(trexp([0 0 1]*0.4), rotz(0.4), 'absTol', 1e-6);
    
    
    %%% SE(3) tests
    
    %% sigma = se(3)
    % pure translation
    tc.verifyEqual(trexp( skewa([1 2 3 0 0 0]) ), transl([1 2 3]), 'absTol', 1e-6);
    tc.verifyEqual(trexp( skewa([0 0 0 0.2 0 0]) ), trotx(0.2), 'absTol', 1e-6);
    tc.verifyEqual(trexp( skewa([0 0 0 0 0.3 0]) ), troty(0.3), 'absTol', 1e-6);
    tc.verifyEqual(trexp( skewa([0 0 0 0 0 0.4]) ), trotz(0.4), 'absTol', 1e-6);
    
    % mixture
    T = transl([1 2 3])*trotx(0.2)*troty(0.3)*trotz(0.4);
    tc.verifyEqual(trexp(logm(T)), T, 'absTol', 1e-6);
    
    %% twist vector
    tc.verifyEqual(trexp( double(Twist(T))), T, 'absTol', 1e-6);
    
    %% (sigma, theta)
    tc.verifyEqual(trexp( skewa([1 0 0 0 0 0]), 2), transl([2 0 0]), 'absTol', 1e-6);
    tc.verifyEqual(trexp( skewa([0 1 0 0 0 0]), 2), transl([0 2 0]), 'absTol', 1e-6);
    tc.verifyEqual(trexp( skewa([0 0 1 0 0 0]), 2), transl([0 0 2]), 'absTol', 1e-6);
    
    tc.verifyEqual(trexp( skewa([0 0 0 1 0 0]), 0.2), trotx(0.2), 'absTol', 1e-6);
    tc.verifyEqual(trexp( skewa([0 0 0 0 1 0]), 0.2), troty(0.2), 'absTol', 1e-6);
    tc.verifyEqual(trexp( skewa([0 0 0 0 0 1]), 0.2), trotz(0.2), 'absTol', 1e-6);
    
    
    %% (twist, theta)
    tc.verifyEqual(trexp(Twist('R', [1 0 0], [0 0 0]).S, 0.3), trotx(0.3), 'absTol', 1e-6);
    
    
    T = transl([1 2 3])*troty(0.3);
    tc.verifyEqual(trexp(logm(T)), T, 'absTol', 1e-6);
    
    tc.verifyError( @() trexp(1), 'SMTB:trexp:badarg')
end

function e2h_test(tc)
    P1 = [1;2; 3];
    P2 = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15];

    tc.verifyEqual( e2h(P1), [P1; 1]);
    tc.verifyEqual( e2h(P2), [P2; ones(1,5)]);
end

function h2e_test(tc)
    P1 = [1;2; 3];
    P2 = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15];

    tc.verifyEqual( h2e(e2h(P1)), P1);
    tc.verifyEqual( h2e(e2h(P2)), P2);
end


function homtrans_test(tc)
    
    P1 = [1;2; 3];
    P2 = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15];
    
    T = eye(4,4);
    tc.verifyEqual( homtrans(T, P1), P1);
    tc.verifyEqual( homtrans(T, P2), P2);
    
    Q = [-2;2; 4];
    T = transl(Q);
    tc.verifyEqual( homtrans(T, P1), P1+Q);
    tc.verifyEqual( homtrans(T, P2), P2+Q);
    
    T = trotx(pi/2);
    tc.verifyEqual( homtrans(T, P1), [P1(1); -P1(3); P1(2)], 'absTol', 1e-6);
    tc.verifyEqual( homtrans(T, P2), [P2(1,:); -P2(3,:); P2(2,:)], 'absTol', 1e-6);
    
    T =  transl(Q)*trotx(pi/2);
    tc.verifyEqual( homtrans(T, P1), [P1(1); -P1(3); P1(2)]+Q, 'absTol', 1e-6);
    tc.verifyEqual( homtrans(T, P2), [P2(1,:); -P2(3,:); P2(2,:)]+Q, 'absTol', 1e-6);

    % projective point case
    P1h = e2h(P1);
    P2h = e2h(P2);
    tc.verifyEqual( homtrans(T, P1h), [P1h(1); -P1h(3); P1h(2); 1]+[Q;0], 'absTol', 1e-6);
    tc.verifyEqual( homtrans(T, P2h), [P2h(1,:); -P2h(3,:); P2h(2,:); ones(1,5)]+[Q;0], 'absTol', 1e-6);

    % sequence case
    TT = cat(3, T, T, T, T, T);
    Q = homtrans(T, TT);
    tc.verifyEqual( Q(:,:,3), T*T);

    % error case
    tc.verifyError( @() homtrans(ones(2,2), P1), 'SMTB:homtrans:badarg')

end

