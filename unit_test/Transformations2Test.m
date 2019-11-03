%% This is for testing the Homogeneous Transformation functions in the robotics Toolbox

function tests = TransformationsTest
    tests = functiontests(localfunctions);
    
    clc
end

function teardownOnce(tc)
    close all
end

%% first of all check we can tell a good matrix from a bad one
function isrot2_test(tc)
    R1 = diag([1 1]);    % proper
    R2 = diag([1 2]);    % not proper
    R3 = diag([2 0.5]);  % not proper
    R4 = diag([1 -1]);  % not proper
    
    % test shapes
    tc.verifyFalse( isrot2(1) )
    tc.verifyFalse( isrot2( zeros(3,3) ) )
    tc.verifyFalse( isrot2( zeros(4,4) ) )
    tc.verifyFalse( isrot2( zeros(2,1) ) )
    tc.verifyFalse( isrot2( zeros(1,2) ) )
    
    % test shapes with validity check
    tc.verifyFalse( isrot2(1, 1) )
    tc.verifyFalse( isrot2( zeros(2,2), 1 ) )
    tc.verifyFalse( isrot2( zeros(4,4) ), 1 )
    tc.verifyFalse( isrot2( zeros(4,1) ), 1 )
    tc.verifyFalse( isrot2( zeros(1,4) ), 1 )
    
    % test 3x3
    tc.verifyTrue( isrot2(R1) )
    tc.verifyTrue( isrot2(R2) )
    tc.verifyTrue( isrot2(R3) )
    tc.verifyTrue( isrot2(R4) )
    
    % test 3x3 with validity check
    tc.verifyTrue( isrot2(R1, 1) )
    tc.verifyFalse( isrot2(R2, 1) )
    tc.verifyFalse( isrot2(R3, 1) )
    tc.verifyFalse( isrot2(R4, 1) )
    
    % vector case
    tc.verifyTrue( isrot2(cat(3, R1, R1, R1)) )
    tc.verifyTrue( isrot2(cat(3, R1, R1, R1), 1) )
    tc.verifyTrue( isrot2(cat(3, R1, R2, R3)) )
    tc.verifyFalse( isrot2(cat(3, R1, R2, R3), 1) )
end

function ishomog2_test(tc)
    T1 = diag([1 1 1]);    % proper
    T2 = diag([1 -1 1]);   % not proper
    T3 = diag([1 2 1]);    % not proper
    T4 = diag([2 0.5 1]);  % not proper
    T5 = diag([1 1 0]);    % not proper
    
    
    % test shapes
    tc.verifyFalse( ishomog2(1) )
    tc.verifyFalse( ishomog2( zeros(2,2) ) )
    tc.verifyFalse( ishomog2( zeros(4,4) ) )
    tc.verifyFalse( ishomog2( zeros(3,1) ) )
    tc.verifyFalse( ishomog2( zeros(1,3) ) )
    
    % test shapes with validity check
    tc.verifyFalse( ishomog2(1, 1) )
    tc.verifyFalse( ishomog2( zeros(2,2), 1 ) )
    tc.verifyFalse( ishomog2( zeros(4,4) ), 1 )
    tc.verifyFalse( ishomog2( zeros(4,1) ), 1 )
    tc.verifyFalse( ishomog2( zeros(1,4) ), 1 )
    
    % test 4x4
    tc.verifyTrue( ishomog2(T1) )
    tc.verifyTrue( ishomog2(T2) )
    tc.verifyTrue( ishomog2(T3) )
    tc.verifyTrue( ishomog2(T4) )
    tc.verifyTrue( ishomog2(T5) )
    
    
    % test 4x4 with validity check
    tc.verifyTrue( ishomog2(T1, 1) )
    tc.verifyFalse( ishomog2(T2, 1) )
    tc.verifyFalse( ishomog2(T3, 1) )
    tc.verifyFalse( ishomog2(T4, 1) )
    tc.verifyFalse( ishomog2(T5, 1) )
    
    
    % vector case
    tc.verifyTrue( ishomog2(cat(3, T1, T1, T1)) )
    tc.verifyTrue( ishomog2(cat(3, T1, T1, T1), 1) )
    tc.verifyTrue( ishomog2(cat(3, T1, T2, T3)) )
    tc.verifyFalse( ishomog2(cat(3, T1, T2, T3), 1) )
end


%% can we convert between rotation matrices and homogeneous coordinate matrices
%    r2t                        - RM to HT
function r2t_test(tc)
    
    % SO(2) case
    R = [1 2; 3 4];
    tc.verifyEqual(r2t(R),...
        [1 2 0; 3 4 0; 0 0 1],'absTol',1e-10);
    
    % sequence case
    Rs = cat(3, R, R, R, R, R);
    Ts = r2t(Rs);
    verifySize(tc, Ts, [3 3 5]);
    tc.verifyEqual(Ts(:,:,2), ...
        [1 2 0; 3 4 0; 0 0 1],'absTol',1e-10);
end


function t2r_test(tc)
    
    % SE(2) case
    T = [1 2 3; 4 5 6; 0 0 1];
    tc.verifyEqual(t2r(T),...
        [1 2; 4 5],'absTol',1e-10);
    
    % sequence case
    Ts = cat(3, T, T, T, T, T);
    Rs = t2r(Ts);
    verifySize(tc, Rs, [2 2 5]);
    tc.verifyEqual(Rs(:,:,2), ...
        [1 2; 4 5],'absTol',1e-10);
end

function rt2tr_test(tc)
    
    R = [1 2; 3 4];
    t = [5; 6];
    
    tc.verifyEqual(rt2tr(R, t),...
        [1 2 5; 3 4 6; 0 0 1],'absTol',1e-10);
    
    % sequence case
    Rs = cat(3, R, 2*R, 3*R);
    ts = cat(2, t, 2*t, 3*t);
    Ts = rt2tr(Rs, ts);
    verifySize(tc, Ts, [3 3 3]);
    tc.verifyEqual(Ts(:,:,1), ...
        [1 2 5; 3 4 6; 0 0 1],'absTol',1e-10);
    tc.verifyEqual(Ts(:,:,2), ...
        [2*[1 2 5; 3 4 6]; 0 0 1],'absTol',1e-10);
    tc.verifyEqual(Ts(:,:,3), ...
        [3*[1 2 5; 3 4 6]; 0 0 1],'absTol',1e-10);
end

function tr2rt_test(tc)
    
    %% SE(2) case
    T = [1 2 3; 4 5 6; 0 0 1];
    [R,t] = tr2rt(T);
    tc.verifyEqual(R, [1 2; 4 5], 'absTol',1e-10);
    tc.verifyEqual(t, [3;6], 'absTol',1e-10);
    
    Ts = cat(3, T, T, T, T, T);
    [Rs,ts] = tr2rt(Ts);
    verifySize(tc, Rs, [2 2 5]);
    verifySize(tc, ts, [5 2]);
    
    tc.verifyEqual(Rs(:,:,2), [1 2; 4 5], 'absTol',1e-10);
    tc.verifyEqual(ts(2,:), [3 6], 'absTol',1e-10);
    
end


%% Constructors

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


%% SO(2)
function rot2_test(tc)
    tc.verifyEqual( rot2(0), eye(2,2), 'absTol', 1e-6);
    tc.verifyEqual( trot2(0), eye(3,3), 'absTol', 1e-6);
    tc.verifyEqual( rot2(0), eye(2,2), 'absTol', 1e-6);
    tc.verifyEqual( trot2(pi/2), [0 -1 0; 1 0 0; 0 0 1], 'absTol', 1e-6);
    
    tc.verifyEqual( trot2(90, 'deg'),[0 -1 0; 1 0 0; 0 0 1], 'absTol', 1e-6);
    tc.verifyEqual( trot2(pi)*trot2(-pi/2), ...
        trot2(pi/2), 'absTol', 1e-6);
    tc.verifyEqual( rot2(pi)*rot2(-pi/2), ...
        rot2(pi/2), 'absTol', 1e-6);
end

%% SE(2)
function SE2_test(tc)
    
    T2 = [1 0 1; 0 1 2; 0 0 1];
    
    % transl(T) -> P
    tc.verifyEqual( transl2(T2), ...
        [1;2], 'absTol',1e-10);
    
    % transl(P) -> T
    tc.verifyEqual( transl2(1, 2), ...
        [1 0 1; 0 1 2; 0 0 1], 'absTol', 1e-6);
    tc.verifyEqual( transl2([2, 3]), ...
        [1 0 2; 0 1 3; 0 0 1], 'absTol', 1e-6);
    
    T3 = transl2([1 2; 3 4; 5 6]);
    verifySize(tc, T3, [3 3 3]);
    tc.verifyEqual(T3(:,:,2), transl2([3 4]))
    
    
    T2f = [1 1 1; 1 1 2; 0 0 1];
    R2 = [1 0 ; 0 1];
    R2f = [1 1 ; 1 1];
    tc.verifyEqual( ishomog2(T2), true);
    tc.verifyEqual( ishomog2(T2,1), true);
    tc.verifyEqual( ishomog2(T2f,1), false);
    tc.verifyEqual( ishomog2(R2), false);
    
    tc.verifyEqual( isrot2(R2), true);
    tc.verifyEqual( isrot2(R2,1), true);
    tc.verifyEqual( isrot2(R2f,1), false);
    tc.verifyEqual( isrot2(T2), false);
    
    
    %     tc.verifyEqual( SE2(2, 3, 0), ...
    %         [1 0 2; 0 1 3; 0 0 1], 'absTol', 1e-6);
    %     tc.verifyEqual( SE2(2, 3, pi/2), ...
    %         transl2(2,3)*trot2(pi/2), 'absTol', 1e-6);
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



function trprint2_test(tc)
    a = transl2([1,2]) * trot2(0.3);
    
    trprint2(a);

    trprint('a') % equivalent to trprint a on the command line
    
    s = evalc( 'trprint2(a)' );
    tc.verifyTrue(isa(s, 'char') );
    tc.verifyEqual(size(s,1), 1);
    
    s = trprint2(a);
    tc.verifyClass(s, 'char');
    tc.verifyEqual(size(s,1), 1);
    
    trprint2(a, 'radian');
    trprint2(a, 'fmt', '%g');
    trprint2(a, 'label', 'bob');
    
    a = cat(3, a, a, a);
    trprint2(a);
    
    s = evalc( 'trprint2(a)' );
    tc.verifyTrue(isa(s, 'char') );
    tc.verifyEqual(size(s,1), 1);
    
    tc.verifyEqual( length(regexp(s, '\n', 'match')), 4);
    
    trprint2(a, 'radian');
    trprint2(a, 'fmt', '%g');
    trprint2(a, 'label', 'bob');
end

function vex_test(tc)
    S = [0 -2; 2 0];
    tc.verifyEqual( vex(S), 2);
    
    S = [0 2; -2 0];
    tc.verifyEqual( vex(S), -2);
end


function skew_test(tc)
    R = skew(3);
    
    tc.verifyTrue( isrot2(R) );  % check size
    
    tc.verifyEqual( norm(R'+ R), 0, 'absTol', 1e-10); % check is skew
    
    tc.verifyEqual( vex(R), 3); % check contents, vex already verified

    tc.verifyError( @() skew([1 2]), 'SMTB:skew:badarg')

end

function vexa_test(tc)
    S = [0 -2 3 ; 2 0 4; 0 0 0];
    tc.verifyEqual( vexa(S), [3 4 2]');
    
    S = [0 2 3 ; -2 0 4; 0 0 0];
    tc.verifyEqual( vexa(S), [3 4 -2]');
end


function skewa_test(tc)
    T = skewa([3 4 5]);
    
    tc.verifyTrue( ishomog2(T) );  % check size
    
    R = T(1:2,1:2);
    tc.verifyEqual( norm(R'+ R), 0, 'absTol', 1e-10); % check is skew
    
    tc.verifyEqual( vexa(T), [3 4 5]'); % check contents, vexa already verified

    tc.verifyError( @() skewa([1 2]), 'SMTB:skewa:badarg')

end


function trexp2_test(tc)
    %unit tests for matrix log stuff
    
    %% so(2)
    tc.verifyEqual( trexp2(skew(0)), rot2(0), 'absTol', 1e-6 );
    tc.verifyEqual( trexp2(skew(pi/2)), rot2(pi/2), 'absTol', 1e-6 );
    tc.verifyEqual( trexp2(skew(-pi/2)), -rot2(pi/2), 'absTol', 1e-6 );
    
    tc.verifyEqual( trexp2(0), rot2(0), 'absTol', 1e-6 );
    tc.verifyEqual( trexp2(pi/2), rot2(pi/2), 'absTol', 1e-6 );
    tc.verifyEqual( trexp2(-pi/2), -rot2(pi/2), 'absTol', 1e-6 );
        
    %% se(2)
    tc.verifyEqual( trexp2(skewa([0 0 0])), transl2([0 0]), 'absTol', 1e-6 );
    tc.verifyEqual( trexp2(skewa([3 0 0])), transl2([3 0]), 'absTol', 1e-6 );
    tc.verifyEqual( trexp2(skewa([0 0 0.3])), trot2(0.3), 'absTol', 1e-6 );
    
    %% errors
    tc.verifyError( @() trexp2(eye(4,4)), 'SMTB:trexp2:badarg');
end

function e2h_test(tc)
    P1 = [1;2];
    P2 = [1 2 3 4 5; 6 7 8 9 10];

    tc.verifyEqual( e2h(P1), [P1; 1]);
    tc.verifyEqual( e2h(P2), [P2; ones(1,5)]);
end

function h2e_test(tc)
    P1 = [1;2];
    P2 = [1 2 3 4 5; 6 7 8 9 10];

    tc.verifyEqual( h2e(e2h(P1)), P1);
    tc.verifyEqual( h2e(e2h(P2)), P2);
end

function homtrans_test(tc)
    
    P1 = [1;2];
    P2 = [1 2 3 4 5; 6 7 8 9 10];
    
    T = eye(3,3);
    tc.verifyEqual( homtrans(T, P1), P1);
    tc.verifyEqual( homtrans(T, P2), P2);
    
    Q = [-2;2];
    T = transl2(Q);
    tc.verifyEqual( homtrans(T, P1), P1+Q);
    tc.verifyEqual( homtrans(T, P2), P2+Q);
    
    T = trot2(pi/2);
    tc.verifyEqual( homtrans(T, P1), [-P1(2); P1(1)], 'absTol', 1e-6);
    tc.verifyEqual( homtrans(T, P2), [-P2(2,:); P2(1,:)], 'absTol', 1e-6);
    
    T =  transl2(Q)*trot2(pi/2);
    tc.verifyEqual( homtrans(T, P1), [-P1(2); P1(1)]+Q, 'absTol', 1e-6);
    tc.verifyEqual( homtrans(T, P2), [-P2(2,:); P2(1,:)]+Q, 'absTol', 1e-6);

    % projective point case
    P1h = e2h(P1);
    P2h = e2h(P2);
    tc.verifyEqual( homtrans(T, P1h), [-P1(2); P1(1); 1]+[Q; 0], 'absTol', 1e-6);
    tc.verifyEqual( homtrans(T, P2h), [-P2(2,:); P2(1,:); ones(1,5)]+[Q; 0], 'absTol', 1e-6);

    % sequence case
    TT = cat(3, T, 2*T, 3*T, 4*T, 5*T);
    Q = homtrans(T, TT);
    tc.verifyEqual( Q(:,:,3), 3*T*T);

    % error case
    tc.verifyError( @() homtrans(7, P1), 'SMTB:homtrans:badarg')
    
end

