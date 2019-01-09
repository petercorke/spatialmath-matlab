% inverse dynamics (recursive Newton-Euler) using spatial vector notation
function  tau = ID( model, q, qd, qdd )
    
    a_grav = SpatialAcceleration([0;0;-9.81;0;0;0]);
    
    for i = 1:model.NB
        [ XJ, S(:,i) ] = jcalc( model.pitch(i), q(i) );
        vJ = SpatialVelocity(S(:,i)*qd(i));
        Xup(i) = XJ*model.Xtree(i);
        if model.parent(i) == 0
            v(i) = vJ;
            a(i) = Xup(i)*(-a_grav) + SpatialAcceleration(S(:,i)*qdd(i));
        else
            v(i) = Xup(i)*v(model.parent(i)) + vJ;
            a(i) = Xup(i)*a(model.parent(i)) ...
                + SpatialAcceleration(S(:,i)*qdd(i)) ...
                + cross(v(i),vJ);
        end
        f(i) = model.I(i)*a(i) + cross( v(i), model.I(i)*v(i) );
    end
    
    for i = model.NB:-1:1
        tau(i,1) = S(:,i)' * double(f(i));
        if model.parent(i) ~= 0
            f(model.parent(i)) = f(model.parent(i)) + Xup(i)*f(i);
        end
    end
end

function  [Xj,S] = jcalc( pitch, q )  %FIXED VW ORDER
    
    % jcalc  Calculate joint transform and motion subspace.
    % [Xj,S]=jcalc(pitch,q) calculates the joint transform and motion subspace
    % matrices for a revolute (pitch==0), prismatic (pitch==inf) or helical
    % (pitch==any other value) joint.  For revolute and helical joints, q is
    % the joint angle.  For prismatic joints, q is the linear displacement.
    
    if pitch == 0				% revolute joint
        Xj = Twist(SE3.Rz(q));
        S = [0;0;0;0;0;1];
    elseif pitch == inf			% prismatic joint
        Xj = Twist(SE3([0 0 q]));
        S = [0;0;1;0;0;0];
    else					% helical joint
        Xj = Twist(SE3.Rz(q) * SE3([0 0 q*pitch]));
        S = [0;0;pitch0;0;1;];
    end
end
