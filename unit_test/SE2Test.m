function tests = SE2Test
    tests = functiontests(localfunctions);
end


% we will assume that the primitives rotx,trotx, etc. all work


function constructor_test(tc)
    
    verifyClass(tc, SE2(), 'SE2');

    %% null
    tc.verifyEqual(SE2().double, eye(3,3));
    
    %% translation only
    t = [1 2];
    tc.verifyEqual(SE2(t).double, transl2(t));
    tc.verifyEqual(SE2(t').double, transl2(t));
    tc.verifyEqual(SE2(t(1), t(2)).double, transl2(t));
    
    %% R
    R = rot2(-pi/2);
    tc.verifyEqual(SE2(R).double, r2t(R));
    
    %% R,t

    tc.verifyEqual(SE2(R,t).double, transl2(t)*r2t(R));
    
    %% T
    T = transl2(1, 2)*trot2(0.3);
    tc.verifyEqual(SE2(T).double, T);
    
    %% x,y,theta
    T = transl2(1, 2)*trot2(0.3);
    tc.verifyEqual(SE2(1, 2, 0.3).double, T);
    tc.verifyEqual(SE2([1, 2, 0.3]).double, T);
    tc.verifyEqual(SE2([1, 2], 0.3).double, T);
        
    %% T
    T = rt2tr(R,t);
    tc.verifyEqual(SE2(T).double, T);
    
    %% copy constructor
    TT = SE2(T);
    tc.verifyEqual(SE2(TT).double, T);
    
    
    %% vectorised versions
    
    T1 = transl2(1,2) * trot2(0.3);
    T2 = transl2(1,-2) * trot2(-0.4);
    
    TT = cat(3, T1, T2, T1, T2);

    tt = SE2(TT);
    tc.verifyEqual(length(tt), size(TT, 3) );
    tc.verifyEqual(tt.T, TT);
    
end

function concat_test(tc)
    x = SE2();
    xx = [x x x x];
    
    tc.verifyClass(xx, 'SE2');
    tc.verifySize(xx, [1 4]);
end

function staticconstructors_test(tc)
    
    %% exponential
    tc.verifyEqual(SE2.exp( skew(0.3) ).R, rot2(0.3), 'AbsTol', 1e-10  );
    
        
    %% exponential
    tc.verifyEqual(SE2.exp(zeros(3,3)).T, eye(3,3), 'AbsTol', 1e-10  );
    t = [1 2];
    tc.verifyEqual(SE2.exp(skewa([t 0])).T, transl2(t), 'AbsTol', 1e-10  );
end



function isa_test(tc)
    
    verifyTrue(tc, SE2.isa(trot2(0)) );
    verifyFalse(tc, SE2.isa(1) )
end

function resulttype_test(tc)
    
    t = SE2();
    verifyClass(tc, t, 'SE2');
    
    verifyClass(tc, t*t, 'SE2');
    
    verifyClass(tc, t/t, 'SE2');
    
    verifyClass(tc, inv(t), 'SE2');
    end

function inverse_test(tc)    
    
    T1 = transl2(1, 2) * trot2(0.3);
    TT1 = SE2(T1);
    
    % test inverse
    tc.verifyEqual(double(TT1.inv()), inv(T1), 'AbsTol', 1e-10  );
    
    tc.verifyEqual(double(TT1*TT1.inv()), eye(3,3), 'AbsTol', 1e-10  );
    tc.verifyEqual(double(TT1.inv()*TT1), eye(3,3), 'AbsTol', 1e-10  );
    
    % vector case
    tc.verifyEqual(double(TT1.inv()*TT1), eye(3,3), 'AbsTol', 1e-10  );
end


function Rt_test(tc)
   
    
    TT1 = SE2.rand
    T1 = TT1.double; R1 = t2r(T1); t1 = transl2(T1);
    
    tc.verifyEqual(TT1.T, T1, 'AbsTol', 1e-10  );
    tc.verifyEqual(TT1.R, R1, 'AbsTol', 1e-10  );
    tc.verifyEqual(TT1.t, t1, 'AbsTol', 1e-10  );
    
    tc.verifyEqual(TT1.transl, t1', 'AbsTol', 1e-10  );
    TT = [TT1 TT1 TT1];
    tc.verifyEqual(TT.transl, [t1 t1 t1]', 'AbsTol', 1e-10  );
end


function arith_test(tc)
    
    R1 = rpy2r( randn(1,3) );  t1 = randn(3,1); T1 = rt2tr(R1, t1);
    R2 = rpy2r( randn(1,3) );  t2 = randn(3,1); T2 = rt2tr(R2, t2);
    
    TT1 = SE2.rand; T1 = TT1.double;
    TT2 = SE2.rand; T2 = TT2.double;

    I = SE2();
    
    
    %% SE2 * SE2 product
    % scalar x scalar
    
    tc.verifyEqual(double(TT1*TT2), T1*T2, 'AbsTol', 1e-10  );
    tc.verifyEqual(double(TT2*TT1), T2*T1, 'AbsTol', 1e-10  );
    tc.verifyEqual(double(TT1*I), T1, 'AbsTol', 1e-10  );
    tc.verifyEqual(double(TT2*I), T2, 'AbsTol', 1e-10  );
    
    % vector x vector
    tc.verifyEqual([TT1 TT1 TT2] * [TT2 TT1 TT1], [TT1*TT2 TT1*TT1 TT2*TT1]);
    
    % scalar x vector
    tc.verifyEqual(TT1 * [TT2 TT1], [TT1*TT2 TT1*TT1]);
    
    % vector x scalar
    tc.verifyEqual([TT1 TT2]*TT2, [TT1*TT2 TT2*TT2]);
    
    %% SE2 * vector product
    vx = [1 0]'; vy = [0 1]'; 

    % scalar x scalar
    
    tc.verifyEqual(TT1*vy, h2e( T1*e2h(vy) ), 'AbsTol', 1e-10);
    
    % vector x vector
    tc.verifyEqual([TT1 TT2] * [vx vy], [h2e(T1*e2h(vx)) h2e(T2*e2h(vy))], 'AbsTol', 1e-10);
    
    % scalar x vector
    tc.verifyEqual(TT1 * [vx vy], h2e( T1*e2h([vx vy]) ), 'AbsTol', 1e-10);
    
    % vector x scalar
    tc.verifyEqual([TT1 TT2 TT1] * vy, [h2e(T1*e2h(vy)) h2e(T2*e2h(vy)) h2e(T1*e2h(vy))], 'AbsTol', 1e-10);
    
end

function function_tests(tc)
    
    % log
    T = SE2.exp([2 3 0.5]);
    tc.verifyEqual(log(T), [0 -0.5 2; 0.5 0 3; 0 0 0], 'AbsTol', 1e-10  );
    
end

function conversions_test(tc)
    
    
    %%  SE2                     convert to SE2 class

    TT = SE2(1, 2, 0.3);
    
    verifyClass(tc, TT.SE3, 'SE3');
    tc.verifyEqual(double(TT.SE3), transl(1, 2, 0) * trotz(0.3), 'AbsTol', 1e-10 );
    
    %% xyt
    tc.verifyEqual(TT.xyt(), [1 2, 0.3]', 'AbsTol', 1e-10);
    
    %% Twist
    T = SE2.exp([2 3 0.5]);
    t = T.Twist();
    verifyInstanceOf(tc, t, 'Twist');
    tc.verifyEqual(t.v, [2 3]', 'AbsTol', 1e-10);
    tc.verifyEqual(t.w, 0.5, 'AbsTol', 1e-10);
    
    %% Lie stuff
    th = 0.3; 
    RR = SO2(th);
    tc.verifyEqual(RR.log, skew(th), 'AbsTol', 1e-10 );

end


function interp_test(tc)
    TT = SE2.rand
    I = SE2;
    
    z = interp(I, TT, 0);
    tc.verifyClass(z, 'SE2')
    
    tc.verifyEqual(double(interp(I, TT, 0)),   double(I), 'AbsTol', 1e-10 );
    tc.verifyEqual(double(interp(I, TT, 1)),   double(TT), 'AbsTol', 1e-4 );
    tc.verifyEqual(double(interp(I, TT, 0.5)), double(trinterp2(TT.T, 0.5)), 'AbsTol', 1e-10  );
    
end



function miscellany_test(tc)
    
    TT = SE2(1, 2, 0.3);
    
    tc.verifyEqual(dim(TT), 3);
        
    tc.verifyEqual(isSE(TT), true );
    
    tc.verifyClass(TT.new, 'SE2');

    tc.verifyClass(SE2.convert(TT), 'SE2');
    tc.verifyClass(SE2.convert(TT.T), 'SE2');
    z = SE2.convert(TT);
    tc.verifyEqual(double(z), double(TT));
    
    z = SE2.convert(TT.T);
    tc.verifyEqual(double(z), TT.T);
    
end


function display_test(tc)
    
    T1 = SE2.rand;
    T2 = SE2.rand
    
    T1.print
    trprint2(T1)   % old style syntax
    
    T1.plot
    
    T1.print
    trprint2(T1)   % old style syntax
    
    T1.plot
    trplot2(T1)   % old style syntax
    
    T1.animate
    T1.animate(T2)
    tranimate2(T1)   % old style syntax
    tranimate2(T1, T2)   % old style syntax
    tranimate2(T1)   % old style syntax
    tranimate2(T1, T2)   % old style syntax
end

