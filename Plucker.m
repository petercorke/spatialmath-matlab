%Plucker Plucker coordinate class
%
% Concrete class to represent a line in Plucker coordinates.
%
% Methods::
% Plucker   Contructor
% -
% pp           Principal point
% ppd          Principal point distance from origin
% point        Generate points on line
% - 
% closest      Closest point on line
% contains     Test if point is on line
% -
% distance     Minimum distance between two lines
% intersects   Intersection point for two lines
% commonperp   Common perpendicular for two lines
% -
% intersect_plane   Intersection points with a plane
% intersect_volume  Intersection points with a volume
% -
% double       Convert to 6-vector
% skew         Convert to 4x4 skew symmetric matrix
% -
% plot         Plot line
% char         Convert to human readable string
% display      Display a human readable string
%
% Operators::
% *       Multiply Plucker matrix by a general matrix
% |       Test if lines are parallel
% ^       Test if lines intersect
% ==      Test if two lines are equivalent
% ~=      Test if lines are not equivalent
%
% Notes::
% - This is reference (handle) class object
% - Link objects can be used in vectors and arrays
%
% References::
% - Ken Shoemake, "Ray Tracing News", Volume 11, Number 1
%   http://www.realtimerendering.com/resources/RTNews/html/rtnv11n1.html#art3
% - Matt Mason lecture notes http://www.cs.cmu.edu/afs/cs/academic/class/16741-s07/www/lectures/lecture9.pdf
%
% Implementation notes::
% - There is a huge variety of notation used across the literature, as well as the ordering
%   of the direction and moment components in the 6-vector


% Copyright (C) 1993-2017, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for MATLAB (RTB).
% 
% RTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% RTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with RTB.  If not, see <http://www.gnu.org/licenses/>.
%
% http://www.petercorke.com

% NOTES
% working: constructor, origin-distance, plane+volume intersect, plot, .L
% method
% TODO
% .L method to skew

classdef Plucker < handle
    
    properties
        w  % direction vector
        v  % moment vector (normal of plane containing line and origin)
    end
    
    properties (Dependent)
        uw  % unit direction vector
    end
    
    methods
        
        
        function pl = Plucker(varargin)
            %Plucker.Plucker Create Plucker line object
            %
            % P = Plucker(L) creates a copy of the Plucker object L
            %
            % P = Plucker(V) creates a Plucker object from V (6x1)
            %
            % P = Plucker(P1, P2) create a Plucker object that represents
            % the line joining the 3D points P1 (3x1) and P2 (3x1).
            %
            % P = Plucker('points', P1, P2) as above.
            %
            % P = Plucker('planes', PI1, PP2) create a Plucker object that represents
            % the line formed by the intersection of two planes PI1, PI2 (4x1).
            %
            % P = Plucker('Pw', P, W) create a Plucker object from a point P (3x1) and
            % direction vector W (3x1).
            %
            % Notes::
            % - Planes are given by the 4-vector [a b c d] to represent ax+by+cz+d=0.

            % simple constructor
            if nargin == 1
                if isvec(varargin{1}, 6)
                    L = varargin{1}; L = L(:);
                    pl.v = L(1:3);
                    pl.w = L(4:6);
                elseif isa(varargin{1}, 'Plucker')
                    L = varargin{1};
                    pl.v = L.v;
                    pl.w = L.w;
                else
                    error('bad arguments to constructor');
                end
                
                return
            end
            
            opt.type = {'points', 'planes', 'Pw', 'wv', 'UQ'};
            [opt,args] = tb_optparse(opt, varargin);
            
            assert(length(args) == 2, 'RTB:Plucker:badarg', 'expecting two vectors');

            % get the two arguments
            A = args{1}; A = A(:);
            B = args{2}; B = B(:);
            

            
            % handle various options
            switch opt.type
                case 'points'
                    assert( isvec(A,3) && isvec(B,3), 'RTB:Plucker:badarg', 'expecting 3-vectors');
                    % compute direction and moment
                    pl.w = A - B;
                    pl.v = cross(A-B, A);
                case 'planes'
                    assert( isvec(A,4) && isvec(B,4), 'RTB:Plucker:badarg', 'expecting 4-vectors');
                    pl.w = cross(A(1:3), B(1:3));
                    pl.v = B(4)*A(1:3) - A(4)*B(1:3);
                case 'wv'
                    assert( isvec(A,3) && isvec(B,3), 'RTB:Plucker:badarg', 'expecting 3-vectors');
                    pl.w = A;
                    pl.v = B;
                case 'Pw'
                    assert( isvec(A,3) && isvec(B,3), 'RTB:Plucker:badarg', 'expecting 3-vectors');
                    %                     pl.P = B;
                    %                     pl.Q = A+B;
                    pl.w = B;
                    pl.v = cross(pl.w, A);
                otherwise
                    error('RTB:Plucker:badarg', 'unknown argument type');
            end
        end
        
        function z = mtimes(a1, a2)
            %Plucker.mtimes Plucker multiplication
            %
            % PL1 * PL2 is the scalar reciprocal product.
            %
            % PL * M is the product of the Plucker skew matrix and M (4xN).
            % M * PL is the product of M (Nx4) and the Plucker skew matrix.
            %
            % Notes::
            % - The * operator is overloaded for convenience.
            % - Multiplication or composition of Plucker lines is not defined.
            %
            % See also Plucker.skew.
            
            if isa(a1, 'Plucker') && isa(a2, 'Plucker')
                % reciprocal product
                z = dot(a1.uw, a2.v) + dot(a2.uw, a1.v);
            elseif isa(a1, 'Plucker') && isfloat(a2)
                assert(numrows(a2) == 4, 'RTB:Plucker:badarg', 'must postmultiply by 4xN matrix');
                z = a1.skew * a2;  % postmultiply by 4xN
            elseif isfloat(a1) && isa(a2, 'Plucker')
                if numcols(a1) == 4
                    z = a1 * a2.skew;  % premultiply by Nx4
                elseif all(size(a1) == [6 6])
                    z = Plucker( a1 * double(a2) ); % premultiply by 6x6 adjoint
                else
                    error('RTB:Plucker:badarg', 'must premultiply by Nx4 matrix');
                end
            end
        end
        
        function x = pp(pl)
            %Plucker.pp Principal point of the line
            %
            % P = PL.pp() is the point on the line that is closest to the origin.
            %
            % Notes::
            % - Same as Plucker.point(0)
            %
            % See also Plucker.ppd, Plucker.point.
            
            x = cross(pl.v, pl.w) / dot(pl.w, pl.w);
        end
        
        function x = double(pl)
            %Plucker.double  Convert Plucker coordinates to real vector
            %
            % PL.double() is a 6x1 vector comprising the moment and direction vectors.
            x = [pl.v; pl.w];
        end
        
        function z = get.uw(pl)
            %Plucker.uw Get the line's direction as a unit vector
            %
            % PL.UW is a unit-vector parallel to the line
            z = unit(pl.w);
        end
        
        function z = skew(pl)
            %Plucker.skew Skew matrix form of the line
            %
            % L = PL.skew() is the Plucker matrix, a 4x4 skew-symmetric matrix
            % representation of the line.
            %
            % Notes::
            % -  For two homogeneous points P and Q on the line, PQ'-QP' is also skew
            %    symmetric.

            
            v = pl.v; w = pl.w;
            
            % the following matrix is at odds with H&Z pg. 72
            z = [
                     0     v(3) -v(2) w(1)
                    -v(3)  0     v(1) w(2)
                     v(2) -v(1)  0    w(3)
                    -w(1) -w(2) -w(3) 0    ];
        end
        
        function z = L(pl)
            warning('RTB:Plucker', 'deprecated: please use skew() method instead');
            z = pl.skew();
        end

        function d = ppd(pl)
            %Plucker.ppd  Distance from principal point to the origin
            %
            % P = PL.ppd() is the distance from the principal point to the origin.
            % This is the smallest distance of any point on the line
            % to the origin.
            %
            % See also Plucker.pp.
            d = sqrt( dot(pl.v, pl.v) / dot(pl.w, pl.w) );
        end
           
        function P = point(L, lambda)
            %Plucker.point Generate point on line
            %
            % P = PL.point(LAMBDA) is a point on the line, where LAMBDA is the parametric
            % distance along the line from the principal point of the line.
            %
            % See also Plucker.pp, Plucker.closest.
            
            P = L.pp + L.uw*lambda(:)';
            %P = bsxfun(@plus, L.P, L.U*lambda(:)');
        end
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  TESTS ON PLUCKER OBJECTS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function t = contains(pl, x)
            %Plucker.contains  Test if point is on the line
            %
            % PL.contains(X) is true if the point X (3x1) lies on the line defined by
            % the Plucker object PL.
            assert( size(x,1) == 3, 'RTB:Plucker: points must have 3 rows');
            t = zeros(1, size(x,2), 'logical');
            for i=1:size(x,2)
                t(i) = norm( cross(x(:,i) - pl.pp, pl.w) ) < 10*eps;
            end
        end
        
        function t = eq(pl1, pl2) 
            %Plucker.eq Test if two lines are equivalent
            %
            % PL1 == PL2 is true if the Plucker objects describe the same line in
            % space.  Note that because of the over parameterization, lines can be
            % equivalent even if they have different dfifferent parameters.
            
            t = abs( 1 - dot(unit(double(pl1)), unit(double(pl2))) ) < 10*eps;
        end
        
        function t = ne(pl1, pl2) 
            %Plucker.ne Test if two lines are not equivalent
            %
            % PL1 ~= PL2 is true if the Plucker objects describe different lines in
            % space.  Note that because of the over parameterization, lines can be
            % equivalent even if they have different dfifferent parameters.
            
            t = abs( 1 - dot(unit(double(pl1)), unit(double(pl2))) ) >= 10*eps;
        end
        
        function v = isparallel(p1, p2)
            %Plucker.isparallel Test if lines are parallel
            %
            % P1.isparallel(P2) is true if the lines represented by Plucker objects P1
            % and P2 are parallel.
            %
            % See also Plucker.or.
            
            v = norm( cross(p1.w, p2.w) ) < 10*eps;
        end
        
        function v = or(p1, p2)
            %Plucker.or Test if lines are parallel
            %
            % P1|P2 is true if the lines represented by Plucker objects P1
            % and P2 are parallel.
            %
            % See also Plucker.isparallel.
            v = isparallel(p1, p2);
        end
        
        function v = mpower(p1, p2)
            %Plucker.mpower Test if lines intersect
            %
            % P1^P2 is true if lines represented by Plucker objects P1
            % and P2 intersect at a point.
            %
            % Notes::
            % - Is false if the lines are equivalent since they would intersect at
            %   an infinite number of points.
            %
            % See also Plucker.parallel.
            v = ~isparallel(p1, p2) && ( abs(p1 * p2) < 10*eps );
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  PLUCKER LINE DISTANCE AND INTERSECTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
       
                
        function p = intersects(p1, p2)
            %Plucker.intersects Find intersection of two lines
            %
            % P = P1.intersects(P2) is the point of intersection (3x1) of the lines
            % represented by Plucker objects P1 and P2.  Is empty matrix if the lines
            % do not intersect, or P1 == P2.
            %
            % Notes::
            % - Is false if the lines are equivalent since they would intersect at
            %   an infinite number of points.
            %
            % See also Plucker.commonperp.

            if p1^p2
                p = -( dot(p1.v,p2.w)*eye(3,3) + p1.w*p2.v' - p2.w*p1.v' ) * unit(cross(p1.w, p2.w));
            else
                p = [];
            end
        end
        
        function l = distance(p1, p2)
            %Plucker.distanve Distance between lines
            %
            % d = P1.distance(P2) is the minimum distance between two lines represented
            % by Plucker objects P1 and P2.
            %
            % Notes::
            % - Works for parallel, skew and intersecting lines.
            if isparallel(p1, p2)
                % lines are parallel
                l = cross(p1.w, p1.v - p2.v * dot(p1.w, p2.w)/ dot(p2.w, p2.w)) / norm(p1.w);
            else
                % lines are not parallel
                if abs(p1 * p2) < 10*eps
                    % lines intersect at a point
                    l = 0;
                else
                    % lines don't intersect, find closest distance
                    l = abs(p1*p2) / norm(cross(p1.w, p2.w))^2;
                end
            end
        end
        
                function [p,dist,lambda] = closest(pl, x)
            %Plucker.closest  Point on line closest to given point
            %
            % P = PL.closest(X) is the coordinate of a point on the line that is
            % closest to the point X (3x1).
            %
            % [P,d] = PL.closest(X) as above but also returns the minimum distance
            % between the point and the line.
            %
            % [P,dist,lambda] = PL.closest(X) as above but also returns the line parameter
            % lambda corresponding to the point on the line.
            %
            % See also Plucker.point.            
            
            % http://www.ahinson.com/algorithms_general/Sections/Geometry/PluckerLine.pdf
            % has different equation for moment, the negative
            
            x = x(:);

            lam = dot(x - pl.pp, pl.uw);
            p = pl.point(lam);  % is the closest point on the line
              
            if nargout > 1
                dist = norm( x - p);
            end
            if nargout > 2
                lambda = lam;
            end
        end
        
        
        function p = commonperp(p1, p2)
            %Plucker.commonperp Common perpendicular to two lines
            %
            % P = P1.commonperp(P2) is a Plucker object representing the common
            % perpendicular line between the lines represented by the Plucker objects
            % P1 and P2.
            %
            % See also Plucker.intersect.
            
            if isparallel(p1, p2)
                % no common perpendicular if lines are parallel
                p = [];
            else
                w = cross(p1.w, p2.w);
                
                v = cross(p1.v, p2.w) - cross(p2.v, p1.w) + ...
                    (p1*p2) * dot(p1.w, p2.w) * unit(cross(p1.w, p2.w));
                
                p = Plucker('wv', w, v);
            end
        end
            

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  PLUCKER LINE DISTANCE AND INTERSECTION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
       
        
        function [p,t] = intersect_plane(L, plane)
            %Plucker.intersect_plane  Line intersection with plane
            %
            % X = PL.intersect_plane(P) is the point where the line intersects the
            % plane P.  Planes are structures with a normal P.n (3x1) and an offset P.p
            % (1x1) such that P.n X + P.p = 0.  X=[] if no intersection.
            %
            % [X,T] = PL.intersect_plane(P) as above but also returns the
            % line parameters (1xN) at the intersection points.
            %
            % See also Plucker.point.
            
            % Line U, V
            % Plane N n
            % (VxN-nU:U.N)
            % Note that this is in homogeneous coordinates.
            %    intersection of plane (n,p) with the line (v,p)
            %    returns point and line parameter
            
            if isstruct(plane)
                N = plane.n;
                n = -dot(plane.n, plane.p);
            else
                N = plane(1:3);
                n = plane(4);
            end
            N = N(:);
            
            den = dot(L.w, N);
            
            if abs(den) > (100*eps)
                %p = -(cross(L.v, N) + n*L.w) / den;
                p = (cross(L.v, N) - n*L.w) / den;
                
                P = L.pp;
                t = dot( P-p, N);
            else
                p = [];
                t = [];
            end
        end

       function [P,lambda] = intersect_volume(line, bounds)
            %PLUCKER.intersect_volume Line intersects plot volume
            %
            % P = PL.intersect_volume(bounds, line) returns a matrix (3xN) with columns
            % that indicate where the line intersects the faces of the plot volume
            % specified in terms of [xmin xmax ymin ymax zmin zmax].  The number of
            % columns N is either 0 (the line is outside the plot volume) or 2 (where
            % the line pierces the bounding volume).  LINE is a structure with elements
            % .p (3x1) a point on the line and .v a vector parallel to the line.
            %
            % [P,T] = PL.intersect_volume(bounds, line) as above but also returns the
            % line parameters (1xN) at the intersection points.
            %
            % See also Plucker.point.
            
            ll = [];
            
            % reshape, top row is minimum, bottom row is maximum
            bounds = reshape(bounds, [2 3]);
            
            for face=1:6
                % for each face of the bounding volume
                %  x=xmin, x=xmax, y=ymin, y=ymax, z=zmin, z=zmax
                
                i = ceil(face/2);  % 1,2,3
                I = eye(3,3);
                plane.n = I(:,i);
                plane.p = [0 0 0]';
                plane.p(i) = bounds(face);
                
                % find where line pierces the plane
                [p,lambda] = line.intersect_plane(plane);

                
                if isempty(p)
                    continue;  % no intersection with this plane
                end
                
%                 fprintf('face %d: n=(%f, %f, %f), p=(%f, %f, %f)\n', face, plane.n, plane.p);
%                 fprintf('      : p=(%f, %f, %f)  ', p)
                
                % find if intersection point is within the cube face
                %  test x,y,z simultaneously
                k = (p' >= bounds(1,:)) & (p' <= bounds(2,:));
                k(i) = [];  % remove the boolean corresponding to current face
                if all(k)
                    % if within bounds, add
                    ll = [ll lambda];
%                     fprintf('  HIT\n');
%                 else
%                     fprintf('\n');
                end
            end
            % put them in ascending order
            ll = sort(ll);
            
            % determine the intersection points from the parameter values
            if isempty(ll)
                P = [];
            else
                P = bsxfun(@plus, line.point(0), line.w*ll);
            end
       end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %  PLOT AND DISPLAY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
        
        function plot(lines, varargin)
            %Plucker.plot Plot a line
            %
            % PL.plot(OPTIONS) plots the Plucker line within the current plot volume.
            %
            % PL.plot(B, OPTIONS) as above but plots within the plot bounds B = [XMIN
            % XMAX YMIN YMAX ZMIN ZMAX].
            %
            % Options::
            % - are passed to plot3.
            %
            % See also plot3.
            bounds = [];
            if nargin > 1
                if all(size(varargin{1}) == [1 6])
                    bounds = varargin{1};
                    varargin = varargin{2:end};
                end
            end
            
            if isempty(bounds)
                bounds = [ get(gca, 'XLim') get(gca, 'YLim') get(gca, 'ZLim')];
            else
                axis(bounds);
            end
            
            %U = pl.Q - pl.P;
            %line.p = pl.P; line.v = unit(U);
            
            ish = ishold();
            hold on
            for pl=lines
                P = pl.intersect_volume(bounds);

                if isempty(P)
                    warning('RTB:Plucker', 'line does not intersect the plot volume');
                else
                    plot3(P(1,:), P(2,:), P(3,:), varargin{:});
                end
            end
            if ~ish
                hold off
            end
        end
        
        
        function display(pl)
            %Plucker.display Display parameters
            %
            % P.display() displays the Plucker parameters in compact single line format.
            %
            % Notes::
            % - This method is invoked implicitly at the command line when the result
            %   of an expression is a Plucker object and the command has no trailing
            %   semicolon.
            %
            % See also Plucker.char.
            loose = strcmp( get(0, 'FormatSpacing'), 'loose');
            if loose
                disp(' ');
            end
            disp([inputname(1), ' = '])
            disp( char(pl) );
        end % display()
        
        function disp(pl)
            disp( char(pl) );
        end
        
        function s = char(pl)
            %Plucker.char Convert to string
            %
            % s = P.char() is a string showing Plucker parameters in a compact single
            % line format.
            %
            % See also Plucker.display.
            s = '';
            for i=1:length(pl)
                
                ps = '{ ';
                ps = [ ps, sprintf('%0.5g  ', pl(i).v) ];
                ps = [ ps(1:end-2), '; '];
                ps = [ ps, sprintf('%0.5g  ', pl(i).w) ];
                ps = [ ps(1:end-2), ' }'];
                if isempty(s)
                    s = ps;
                else
                    s = char(s, ps);
                end
            end
        end
        
        
%         function z = side(pl1, pl2)
%             %Plucker.side Plucker side operator
%             %
%             % X = SIDE(P1, P2) is the side operator which is zero whenever
%             % the lines P1 and P2 intersect or are parallel.
%             %
%             % See also Plucker.or.
%             
%             if ~isa(pl2, 'Plucker')
%                 error('RTB:Plucker:badarg', 'both arguments to | must be Plucker objects');
%             end
%             L1 = pl1.line(); L2 = pl2.line();
%             
%             z = L1([1 5 2 6 3 4]) * L2([5 1 6 2 4 3])';
%         end

%         
%         function z = intersect(pl1, pl2)
%             Plucker.intersect  Line intersection
%             
%             PL1.intersect(PL2) is zero if the lines intersect.  It is positive if PL2
%             passes counterclockwise and negative if PL2 passes clockwise.  Defined as
%             looking in direction of PL1
%             
%                                        ---------->
%                            o                o
%                       ---------->
%                      counterclockwise    clockwise
%             
%             z = dot(pl1.w, pl1.v) + dot(pl2.w, pl2.v);
%         end
        
    end % methods
end % class


