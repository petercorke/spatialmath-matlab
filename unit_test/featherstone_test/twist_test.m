a = Twist(SE3([1 2 3]))
%a.SE
% a.Ad
x = inv(a.SE); 
x.Ad

Xtrans([1 2 3])

a = Twist(SE3.Rx(pi/2))
%a.SE
x = inv(a.SE); 
x.Ad

Xrotx(pi/2)