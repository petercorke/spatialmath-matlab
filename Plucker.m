%Plucker Plucker coordinate class
%
% Concrete class to represent a line in Plucker coordinates.
%
% Methods::
% line    Return Plucker line coordinates (1x6)
% side    Side operator
%
% Operators::
% *       Multiple Plucker matrix by a general matrix
% |       Side operator
%
%
% Notes::
% - This is reference class object
% - Link objects can be used in vectors and arrays
%
% References::
% - Ken Shoemake, "Ray Tracing News", Volume 11, Number 1

% Copyright (C) 1993-2014, by Peter I. Corke
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
        P  % a point on the line (if given)
        Q  % 
        U  % direction vector  also w
        V  % moment vector
    end
    
%     properties(Dependent=true)
%         v
%         w
%     end
    
    methods
        
        function w = w(pl)
            w = pl.U;
        end
        
                function v = v(pl)
            v = pl.V;
        end
        
        function pl = Plucker(varargin)
            %Plucker.Plucker Create Plucker object
            %
            % P = Plucker(P1, P2) create a Plucker object that represents
            % the line joining the 3D points P1 (3x1) and P2 (3x1).
            
            opt.type = {'points', 'planes', 'Pw', 'UV', 'UQ'};
            [opt,args] = tb_optparse(opt, varargin);
            
            if length(args) ~= 2
                error('RTB:Plucker:badarg', 'expecting two vectors');
            end
            
            A = args{1}; A = A(:);
            B = args{2}; B = B(:);

   
            switch opt.type
                    
                case 'points'
                    if ~isvec(A) || ~isvec(B)
                        error('RTB:Plucker:badarg', 'expecting 3-vectors');
                    end
                    % stash the points
                    
                    pl.P = B;
                    pl.Q = A;
                    % compute direction and moment
                    pl.U = B - A;
                    pl.V = cross(B, A);
                case 'planes'
                    if ~isvec(A,4) || ~isvec(B,4)
                        error('RTB:Plucker:badarg', 'expecting 4-vectors');
                    end
                    pl.U = cross(A(1:3), B(1:3));
                    pl.V = B(4)*A(1:3) - A(4)*B(1:3);
                case 'UV'
                    if ~isvec(A) || ~isvec(B)
                        error('RTB:Plucker:badarg', 'expecting 3-vectors');
                    end
                    pl.U = A;
                    pl.V = B;
               case 'Pw'
                    if ~isvec(A) || ~isvec(B)
                        error('RTB:Plucker:badarg', 'expecting 3-vectors');
                    end
                    pl.P = B;
                    pl.P = A+B;
                    pl.U = B;
                    pl.V = cross(B, A);
                case 'UQ'
                    if ~isvec(A) || ~isvec(B)
                        error('RTB:Plucker:badarg', 'expecting 3-vectors');
                    end
                    pl.Q = B;
                    pl.P = A+B;
                    pl.U = A;
                    pl.V = cross(A, B);
                otherwise
            end
            

        end
        
        function p = origin_closest(pl)
            p = cross(pl.V, pl.U) / dot(pl.U, pl.U);
        end
        
        function d = origin_distance(pl)
            d = sqrt( dot(pl.V, pl.V) / dot(pl.U, pl.U) );
        end
        
        function P = point(L, lambda)
            
            
           P = bsxfun(@plus, L.P, L.U*lambda(:)');
        end

        
        function z = L(pl)
            % create the Plucker matrix
            Ph = e2h(pl.P);
            Qh = e2h(pl.Q);
            z = Ph*Qh' - Qh*Ph';
        end
        
        function v = line(pl)
            %Plucker.double Plucker line coordinates
            %
            % P.line() is a 6-vector representation of the Plucker
            % coordinates of the line.
            
%             L = pl.L;
%             v = [L(2,1) L(3,1) L(4,1) L(4,3) L(2,4) L(3,2)];
             v = [pl.U; pl.V]';
        end
        
        function z = mtimes(a1, a2)
            %Plucker.mtimes Plucker composition
            %
            % P * M is the product of the Plucker matrix and M (4xN).
            %
            % M * P is the product of M (Nx4) and the Plucker matrix.
            
            if isa(a1, 'Plucker')
                if numrows(a2) ~= 4
                    error('RTB:Plucker:badarg', 'must postmultiply by 4xN matrix');
                end
                z = a1.L * a2;
            elseif isa(a2, 'Plucker')
                if numcols(a1) ~= 4
                    error('RTB:Plucker:badarg', 'must premultiply by Nx4 matrix');
                end
                z = a1 * a2.L;
            end
        end
        
        function z = or(pl1, pl2)
            % P1 | P2 is the side operator which is zero whenever
            % the lines P1 and P2 intersect or are parallel.
            z = side(pl1, pl2);
        end
        
        function z = side(pl1, pl2)
            %Plucker.mtimes Side operator
            %
            % SIDE(P1, P2) is the side operator which is zero whenever
            % the lines P1 and P2 intersect or are parallel.
            if ~isa(pl2, 'Plucker')
                error('RTB:Plucker:badarg', 'both arguments to | must be Plucker objects');
            end
            L1 = pl1.line(); L2 = pl2.line();
            
            z = L1([1 5 2 6 3 4]) * L2([5 1 6 2 4 3])';
        end
        
        function z = intersect(pl1, pl2)
            % is zero if they intersect > 0 if pass counterclockwise
            z = dot(pl1.w, pl1.v) 
            
%                      ---------->
%          o                o
%     ---------->
%   counterclockwise    clockwise
  + dot(pl2.w, pl2.v);
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
            p = (cross(pl.V, P(1:3)) - P(4)*pl.U) / dot(pl.U, P(1:3));
        end
        
                
        function [p,t] = intersect_plane(L, plane)
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
            
            den = dot(L.U, N);
            
           if abs(den) > (100*eps)
                p = (cross(L.v, N) - n*L.w) / den;
                t = dot( p-L.P, N);
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
        
        function plot(pl, varargin)
            
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
            
                U = pl.Q - pl.P;
            %line.p = pl.P; line.v = unit(U);
            P = pl.intersect_volume(bounds)

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
                ps = [ ps, sprintf('%0.5g  ', pl.V) ];
                ps = [ ps(1:end-2), '; '];
                ps = [ ps, sprintf('%0.5g  ', pl.U) ];
                ps = [ ps(1:end-2), ' }'];
                if isempty(s)
                    s = ps;
                else
                    s = char(s, ps);
                end
            end
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
                [p,t] = intersect_plane(line, plane);
                
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
            P = bsxfun(@plus, line.P, line.U*tt);
        end
        
    end % methods
end % class

 
        