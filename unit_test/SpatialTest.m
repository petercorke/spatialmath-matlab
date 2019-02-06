%% This is for testing the Trajectory Generation functions in the robotics Toolbox
function tests = SpatialTest
    tests = functiontests(localfunctions);
    clc
end


function velocity_test(tc)
    a = SpatialVelocity([1:6]');
    tc.verifyClass(a, 'SpatialVelocity')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifyTrue( isa(a, 'SpatialM6') )
    tc.verifySize(a, [1 1]);
    
    a = SpatialVelocity([1:6]);
    tc.verifyClass(a, 'SpatialVelocity')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifySize(a, [1 1]);
    
    s = char(a);
    tc.verifyTrue(ischar(s));
    tc.verifyTrue(size(s,1) == 1);
    tc.verifyTrue(startsWith(s, 'SpatialVelocity'))
    
    v = rand(6,1);
    a = SpatialVelocity(v);
    tc.verifyEqual( double(a), v, 'AbsTol', 1e-10);
    
    r = rand(6,10);
    a = SpatialVelocity(r);
    tc.verifyClass(a, 'SpatialVelocity')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifySize(a, [1 10]);
    
    tc.verifyEqual( double(a), r, 'AbsTol', 1e-10);
    
    s = char(a);
    tc.verifyTrue(ischar(s));
    tc.verifyTrue(size(s,1) == 10);
    
    
end


function vacceleration_test(tc)
    a = SpatialAcceleration([1:6]');
    tc.verifyClass(a, 'SpatialAcceleration')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifyTrue( isa(a, 'SpatialM6') )
    tc.verifySize(a, [1 1]);
    
    a = SpatialAcceleration([1:6]);
    tc.verifyClass(a, 'SpatialAcceleration')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifySize(a, [1 1]);
    
    s = char(a);
    tc.verifyTrue(ischar(s));
    tc.verifyTrue(size(s,1) == 1);
    tc.verifyTrue(startsWith(s, 'SpatialAcceleration'))
    
    v = rand(6,1);
    a = SpatialAcceleration(v);
    tc.verifyEqual( double(a), v, 'AbsTol', 1e-10);
    
    r = rand(6,10);
    a = SpatialAcceleration(r);
    tc.verifyClass(a, 'SpatialAcceleration')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifySize(a, [1 10]);
    
    tc.verifyEqual( double(a), r, 'AbsTol', 1e-10);
    
    s = char(a);
    tc.verifyTrue(ischar(s));
    tc.verifyTrue(size(s,1) == 10);
    
    
end

function force_test(tc)
    
    a = SpatialForce([1:6]');
    tc.verifyClass(a, 'SpatialForce')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifyTrue( isa(a, 'SpatialF6') )
    tc.verifySize(a, [1 1]);
    
    a = SpatialForce([1:6]);
    tc.verifyClass(a, 'SpatialForce')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifyTrue( isa(a, 'SpatialF6') )
    tc.verifySize(a, [1 1]);
    
    s = char(a);
    tc.verifyTrue(ischar(s));
    tc.verifyTrue(size(s,1) == 1);
    tc.verifyTrue(startsWith(s, 'SpatialForce'))
    
    v = rand(6,1);
    a = SpatialForce(v);
    tc.verifyEqual( double(a), v, 'AbsTol', 1e-10);
    
    r = rand(6,10);
    a = SpatialForce(r);
    tc.verifyClass(a, 'SpatialForce')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifyTrue( isa(a, 'SpatialF6') )
    tc.verifySize(a, [1 10]);
    
    tc.verifyEqual( double(a), r, 'AbsTol', 1e-10);
    
    s = char(a);
    tc.verifyTrue(ischar(s));
    tc.verifyTrue(size(s,1) == 10);
    
end

function momentum_test(tc)
    
    a = SpatialMomentum([1:6]');
    tc.verifyClass(a, 'SpatialMomentum')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifyTrue( isa(a, 'SpatialF6') )
    tc.verifySize(a, [1 1]);
    
    a = SpatialMomentum([1:6]);
    tc.verifyClass(a, 'SpatialMomentum')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifyTrue( isa(a, 'SpatialF6') )
    tc.verifySize(a, [1 1]);
    
    s = char(a);
    tc.verifyTrue(ischar(s));
    tc.verifyTrue(size(s,1) == 1);
    tc.verifyTrue(startsWith(s, 'SpatialMomentum'))
    
    v = rand(6,1);
    a = SpatialMomentum(v);
    tc.verifyEqual( double(a), v, 'AbsTol', 1e-10);
    
    r = rand(6,10);
    a = SpatialMomentum(r);
    tc.verifyClass(a, 'SpatialMomentum')
    tc.verifyTrue( isa(a, 'SpatialVec6') )
    tc.verifyTrue( isa(a, 'SpatialF6') )
    tc.verifySize(a, [1 10]);
    
    tc.verifyEqual( double(a), r, 'AbsTol', 1e-10);
    
    s = char(a);
    tc.verifyTrue(ischar(s));
    tc.verifyTrue(size(s,1) == 10);
    
end

function arith_test(tc)
    % just test SpatialVelocity since all types derive from same superclass
    
    r1 = rand(6,1); r2 = rand(6,1);
    a1 = SpatialVelocity(r1);
    a2 = SpatialVelocity(r2);
    
    tc.verifyEqual( double(a1+a2), r1+r2, 'AbsTol', 1e-10);
    tc.verifyEqual( double(a1-a2), r1-r2, 'AbsTol', 1e-10);
    tc.verifyEqual( double(-a1), -r1, 'AbsTol', 1e-10);
end
