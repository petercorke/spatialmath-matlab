%Plucker Plucker coordinate class
%
% Concrete class to represent a line in Plucker coordinates.
%
% Methods::
% line    Return Plucker line coordinates (1x6)
% side    Side operator
% origin_closest
% origin_distance
% distance
% mindist
% point
% pp
% L
% intersect
% Operators::
% *       Multiply Plucker matrix by a general matrix
% |       Side operator
%
%
% Notes::
% - This is reference class object
% - Link objects can be used in vectors and arrays
%
% References::
% - Ken Shoemake, "Ray Tracing News", Volume 11, Number 1
%   http://www.realtimerendering.com/resources/RTNews/html/rtnv11n1.html#art3

% Implementation notes:
% lots of different notation about
%  Shoemaker U:V is (w,v)


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
        v  % moment vector
    end
    
    
    methods
        
        
        function pl = Plucker(varargin)
            %Plucker.Plucker Create Plucker object
            %
            % P = Plucker(P1, P2) create a Plucker object that represents
            % the line joining the 3D points P1 (3x1) and P2 (3x1).
            %
            % P = Plucker('points', P1, P2) as above.
            %
            % P = Plucker('planes', PL1, PL2) create a Plucker object that represents
            % the line formed by the intersection of two planes PL1, PL2 (4x1).
            %
            % P = Plucker('wv', W, V) create a Plucker object from its direction W (3x1) and
            % moment vectors V (3x1).
            %
            % P = Plucker('Pw', P, W) create a Plucker object from a point P (3x1) and
            % direction vector W (3x1).

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
                    pl.v = cross(pl.w, A);
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
                    pl.v = cross(B, A);
                otherwise
                    error('RTB:Plucker:badarg', 'unknown argument type');
            end
        end
        
        function p = origin_closest(pl)
            %Plucker.origin_closest  Point on line closest to the origin
            %
            % P = PL.origin_closest() is the coordinate of a point on the line that is
            % closest to the origin.
            %
            % See also Plucker.origin_distance.
            
            p = cross(pl.v, pl.w) / dot(pl.w, pl.w);
        end
        
        function d = origin_distance(pl)
            %Plucker.origin_distance  Smallest distance from line to the origin
            %
            % P = PL.origin_distance() is the smallest distance of a point on the line
            % to the origin.
            %
            % See also Plucker.origin_closest.
            d = sqrt( dot(pl.v, pl.v) / dot(pl.w, pl.w) );
        end
        
        function [p,dist] = closest(pl, x)
            %Plucker.closest  Point on line closest to given point
            %
            % P = PL.closest(X) is the coordinate of a point on the line that is
            % closest to the point X (3x1).
            %
            % [P,D] = PL.closest(X) as above but also returns the closest distance.
            %
            % See also Plucker.origin_closest.            
            
            % http://www.ahinson.com/algorithms_general/Sections/Geometry/PluckerLine.pdf
            % has different equation for moment, the negative
            w = unit(pl.w);
            pp = pl.pp;
            
            d = dot(x.w, w);
            p = pp + d*w;
            
            if nargout == 2
                dist = norm(p - x.w);
            end
        end
        
        function d = mindist(pl1, pl2)
            %Plucker.mindist  Minimum distance between two lines
            %
            % D = PL1.mindist(PL2) is the minimum distance between two Plucker lines
            % PL1 and PL2.
            
            d = dot(pl1.w, pl2.v) + dot(pl2.w, pl1.v);
        end
        
        function P = point(L, lambda)
            %Plucker.point Point on line
            %
            % P = PL.point(L) is a point on the line, where L is the parametric
            % distance along the line from the principal point of the line.
            %
            % See also Plucker.pp.
            
            P = cross(L.w, L.v) / dot(L.w, L.w) + L.w*lambda(:)';
            %P = bsxfun(@plus, L.P, L.U*lambda(:)');
        end
        
        function x = pp(pl)
            %Plucker.pp Principal point of the line
            %
            % P = PL.pp() is a point on the line.
            %
            % Notes::
            % - Same as Plucker.point(0)
            %
            % See also Plucker.point.
            
            % principal point, for lambda = 0
            x = cross( pl.v, unit(pl.w) );
        end
        
        function z = L(pl)
            %Plucker.L Skew matrix form of the line
            %
            % L = PL.L() is the Plucker matrix, a 4x4 skew-symmetric matrix
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
        
        function v = line(pl)
            %Plucker.double Plucker line coordinates
            %
            % P.line() is a 6-vector representation of the Plucker
            % coordinates of the line.
            %
            % See also Plucker.v, Plucker.w.
            
            %             L = pl.L;
            %             v = [L(2,1) L(3,1) L(4,1) L(4,3) L(2,4) L(3,2)];
            v = [pl.w; pl.v]';
        end
        
        
        
        function z = mtimes(a1, a2)
            %Plucker.mtimes Plucker composition
            %
            % PL * M is the product of the Plucker matrix and M (4xN).
            %
            % M * PL is the product of M (Nx4) and the Plucker matrix.
            
            if isa(a1, 'Plucker')
                assert(numrows(a2) == 4, 'RTB:Plucker:badarg', 'must postmultiply by 4xN matrix');
                z = a1.L * a2;
            elseif isa(a2, 'Plucker')
                assert(numrows(a2) == 4, 'RTB:Plucker:badarg', 'must premultiply by Nx4 matrix');
                z = a1 * a2.L;
            end
        end
        
        function z = or(pl1, pl2)
            %Plucker.or Operator form of side operator
            %
            % P1 | P2 is the side operator which is zero whenever
            % the lines P1 and P2 intersect or are parallel.
            %
            % See also Plucker.side.
            
            z = side(pl1, pl2);
        end
        
        function z = side(pl1, pl2)
            %Plucker.side Plucker side operator
            %
            % X = SIDE(P1, P2) is the side operator which is zero whenever
            % the lines P1 and P2 intersect or are parallel.
            %
            % See also Plucker.or.
            
            if ~isa(pl2, 'Plucker')
                error('RTB:Plucker:badarg', 'both arguments to | must be Plucker objects');
            end
            L1 = pl1.line(); L2 = pl2.line();
            
            z = L1([1 5 2 6 3 4]) * L2([5 1 6 2 4 3])';
        end
        
        function z = intersect(pl1, pl2)
            %Plucker.intersect  Line intersection
            %
            % PL1.intersect(PL2) is zero if the lines intersect.  It is positive if PL2
            % passes counterclockwise and negative if PL2 passes clockwise.  Defined as
            % looking in direction of PL1
            %
            %                            ---------->
            %                o                o
            %           ---------->
            %          counterclockwise    clockwise
            
            z = dot(pl1.w, pl1.v) + dot(pl2.w, pl2.v);
        end
        
        % operator tests for
        % identical ==
        % parallel ||
        % skewed - find closest point
        % intersecting - find
        % distance from origin  d^2 = (V.V)/U.U
        % closest pnt to origin (VxU:U.U)
        
        % planes E.P + e = 0, F.P +f = 0
        % L = ExF : fE - eF
        
        % line plane intersection
        %(VxN-nU:U.N)
        
        % common plane
        % (UxN
        
        function p = plane_intersect(pl, P)
            P = P(:);
            p = (cross(pl.v, P(1:3)) - P(4)*pl.w) / dot(pl.w, P(1:3));
        end
        
        
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
            N= N(:);
            
            den = dot(L.w, N);
            
            if abs(den) > (100*eps)
                %p = -(cross(L.v, N) + n*L.w) / den;
                p = (-cross(L.v, N) - n*L.w) / den;
                
                P = L.point(0);
                t = dot( p-P, N);
            else
                p = [];
                t = [];
            end
            
            
            %             if abs(dvn) > (100*eps)
            %                 t = dot( p-L.P, n) / dot(L.v, n);
            %                 p = L.P + t*L.v;
            %             else
            %                 p = [];
            %                 t = [];
            %             end
        end

       function [P,t] = intersect_volume(line, bounds)
            %PLUCKER.intersect_volume Line intersects plot volume
            %
            % P = PL.intersect_volume(bounds, line) returns a matrix (3xN) with columns that
            % indicate where the line intersects the faces of the plot volume specified
            % in terms of [xmin xmax ymin ymax zmin zmax].  The number of columns N is
            % either 0 (the line is outside the plot volume) or 2.  LINE is a structure
            % with elements .p (3x1) a point on the line and .v a vector parallel to
            % the line.
            %
            % [P,T] = PL.intersect_volume(bounds, line) as above but also returns the
            % line parameters (1xN) at the intersection points.
            %
            % See also Plucker.point.
            
            tt = [];
            
            % reshape, top row is minimum, bottom row is maximum
            bounds = reshape(bounds, [2 3]);
            
            for face=1:6
                % for each face of the bounding volume
                
                i = ceil(face/2);  % 1,2,3
                I = eye(3,3);
                plane.n = I(:,i);
                plane.p = [0 0 0]';
                plane.p(i) = bounds(face);
                [p,t] = line.intersect_plane(plane);
                
                if isempty(p)
                    continue;  % no intersection
                end
                
                k = (p' > bounds(1,:)) & (p' < bounds(2,:));
                k(i) = [];
                if all(k)
                    tt = [tt t];
                end
            end
            % put them in ascending order
            tt = sort(tt);
            
            % determine the intersection points from the parameter values
            P = bsxfun(@plus, line.point(0), line.w*tt);
            
                end
        
        function x = double(pl)
            %Plucker.double  Convert Plucker coordinates to real vector
            %
            % PL.double() is a 6x1 vector comprising the moment and direction vectors.
            x = [pl.v; pl.w];
        end
        
        function plot(pl, varargin)
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
            P = pl.intersect_volume(bounds);
            
            plot3(P(1,:), P(2,:), P(3,:), varargin{:})
            
            %             plot3(P(1,1), P(2,1), P(3,1), 'bx')
            %             plot3(P(1,2), P(2,2), P(3,2), 'bx')
            %             %U = U/2; V = V/2;
            %             L = [cross(V,U); dot(U,U)]
            %             LL = L / L(4);
            %             plot3(LL(1), LL(2), LL(3), 'o');
            %             plot3([0 LL(1)], [0 LL(2)], [0 LL(3)], 'r');
        end
        
        %         @ L = {U:V}, with 3-tuples U and V, with U.V = 0, and with U non-null.
        %   @ L = {P-Q:PxQ}, for P and Q distinct points on L, and line is directed Q->P.
        %   @ L = {U:UxQ}, for U the direction of L and Q a point on L.
        %   @ L = {qP-pQ:PxQ}, for (P:p) and (Q:q) distinct homogeneous points on L.
        %   @ L = {ExF:fE-eF}, for [E:e] and [F:f] distinct planes containing L.
        %   @ {U1:V1} =? s{U2:V2} tests if L1 = {U1:V1} equals L2 = {U2:V2}.
        %   @ s > 0 if L1 and L2 have same orientation.
        %   @ (V.V)/(U.U) is the minimum squared distance of L from the origin.
        %   @ (VxU:U.U) is the point of L closest to the origin.
        %
        %    @ (VxN-Un:U.N) is the point where L intersects plane [N:n] not parallel to L.
        %   @ [UxP-Vw:V.P] is the plane containing L and point (P:w) not on L.
        
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
                ps = [ ps, sprintf('%0.5g  ', pl.v) ];
                ps = [ ps(1:end-2), '; '];
                ps = [ ps, sprintf('%0.5g  ', pl.w) ];
                ps = [ ps(1:end-2), ' }'];
                if isempty(s)
                    s = ps;
                else
                    s = char(s, ps);
                end
            end
        end
        

        
    end % methods
end % class


