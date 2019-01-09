%% create robot model
function tau = idyn_roy(q, qd, qdd)
% copied from planarN
n = 2;
robot.NB = n;
robot.parent = [0:n-1];
robot.pitch = zeros(1,n);  % revolute joints


robot.Xtree(1) = Twist(inv(SE3.Rx(pi/2)));  % joint 1 at origin, z in horizontal plane
robot.Xtree(2) = Twist(inv(SE3([1 0 0])));

robot.I(1) = SpatialInertia( 1, [0.5 0 0], diag([0 0 0]) );
robot.I(2) = SpatialInertia( 1, [0.5 0 0], diag([0 0 0]) );

%% compute inverse dynamics
tau = ID( robot, q, qd, qdd)' %, f_ext, grav_accn )
end