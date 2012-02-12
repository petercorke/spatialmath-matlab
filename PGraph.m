%PGraph Simple graph class
%
%   g = PGraph()        create a 2D, planar, undirected graph
%   g = PGraph(n)       create an n-d, undirected graph
%
% Graphs::
%   - are undirected
%   - are symmetric cost edges (A to B is same cost as B to A)
%   - are embedded in coordinate system
%   - have no loops (edges from A to A)
%   - vertices are represented by integer ids, vid
%   - edges are represented by integer ids, eid
%
%   Graph connectivity is maintained by a labeling algorithm and this
%   is updated every time an edge is added.
%
% Methods::
%
% Constructing the graph::
%   g.add_node(coord)      add vertex, return vid
%   g.add_node(coord, v)   add vertex and edge to v, return vid
%   g.add_edge(v1, v2)     add edge from v1 to v2, return eid
%
%   g.clear()              remove all nodes and edges from the graph
%
% Information from graph::
%   g.edges(e)             return vid for edge
%   g.cost(e)              return cost for edge list
%   g.coord(v)             return coordinate of node v
%   g.neighbours(v)        return vid for edge
%   g.component(v)         return component id for vertex
%   g.connectivity()       return number of edges for all nodes
%
%   g.plot()               set goal vertex for path planning
%   g.pick()               return vertex id closest to picked point
%
%   char(g)                display summary info about the graph
%
% Planning paths through the graph::
%   g.goal(v)              set goal vertex, and plan paths
%   g.next(v)              return d of neighbour of v closest to goal
%   g.path(v)              return list of nodes from v to goal
%
% Graph and world points::
%   g.distance(v1, v2)     distance between v1 and v2 as the crow flies
%   g.closest(coord)       return vertex closest to coord
%   g.distances(coord)     return sorted distances from coord and vertices
%
% To change the distance metric create a subclass of PGraph and override the 
% method distance_metric().
%
% Object properties (read/write)::
%   g.n            number of nodes
%

% Peter Corke 8/2009.

% TODO: add support for arbitrary distance function
%       add A* for distance between 2 vertices
%       be able to delete nodes, must update connectivity
%       add adjacency matrix

classdef PGraph < handle

    properties
        vertexlist        % vertex coordinates, columnwise, vertex number is the column number
        edgelist        % 2xNe matrix, each column is vertex index of edge start and end
        edgelen         % length (cost) of this edge    

        curLabel        % current label
        ncomponents     % number of components
        labels          % label of each vertex
        labelset        % set of labels

        goaldist        % distance from goal, after planning
        
        n               % number of nodes Nv
        
        userdata        % per vertex data, cell array
        ndims           % number of coordinate dimensions, height of vertices matrix
        verbose
        measure         % distance measure: 'Euclidean', 'SE2'
    end

    methods

        function g = PGraph(ndims, varargin)
        %PGraph.PGraph Graph class constructor
        %
        % G = PGraph(D, OPTIONS) returns a graph object embedded 
        % in D dimensions.
        %
        % Options::
        %  'distance', M   Use the distance metric M for path planning
        %  'verbose'       Specify verbose operation
        %
        % Note::
        % - The distance metric is either 'Euclidean' or 'SE2' which is the sum of
        %   the squares of the difference in position and angle modulo 2pi.

            if nargin < 1
                ndims = 2;  % planar by default
            end
            g.ndims = ndims;
            opt.distance = 'Euclidean';
            opt = tb_optparse(opt, varargin);
            
            g.clear();
            g.verbose = opt.verbose;
            g.measure = opt.distance;
        end

        function clear(g)
        %PGraph.CLEAR Clear the graph
        %
        % G.CLEAR() removes all nodes and edges.
            g.labelset = zeros(1, 0);
            g.labels = zeros(1, 0);
            g.edgelist = zeros(2, 0);
            g.edgelen = zeros(1, 0);
            g.vertexlist = zeros(g.ndims, 0);

            g.ncomponents = 0;
            g.curLabel = 0;
            g.n = 0;
        end

        function v = add_node(g, coord, varargin)
        %PGraph.add_node Add a node to the graph
        %
        % V = G.add_node(X) adds a node with coordinate X, where X is Dx1, and
        % returns the node id V.
        %
        % V = G.add_node(X, V) adds a node with coordinate X and connected to
        % node V by an edge.
        %
        % V = G.add_node(X, V, C) adds a node with coordinate X and connected to
        % node V by an edge with cost C.

            if length(coord) ~= g.ndims
                error('coordinate length different to graph coordinate dimensions');
            end
            
            % append the coordinate as a column in the vertex matrix
            g.vertexlist = [g.vertexlist coord(:)];
            v = numcols(g.vertexlist);
            g.labels(v) = g.newlabel();
            if nargin > 2
                g.add_edge(v, varargin{:});
            end
            g.n = g.n + 1;
            if g.verbose
                fprintf('add node (%d) = ', v);
                fprintf('%f ', coord);
                fprintf('\n');
            end
        end
        
        function add_edge(g, v1, v2, d)
        %PGraph.add_edge Add an edge to the graph
        %
        % E = G.add_edge(V1, V2) add an edge between nodes with id V1 and V2, and
        % returns the edge id E.
        %
        % E = G.add_edge(V1, V2, C) add an edge between nodes with id V1 and V2 with
        % cost C.
            if g.verbose
                fprintf('add edge %d -> %d\n', v1, v2);
            end
            for vv=v2(:)'
                g.edgelist = [g.edgelist [v1; vv]];
                if nargin < 4
                    d = g.distance(v1, vv);
                end
                g.edgelen = [g.edgelen d];
                if g.labels(vv) ~= g.labels(v1)
                    g.merge(g.labels(vv), g.labels(v1));
                end
            end
        end

        function c = component(g, v)
            c = [];
            for vv=v
                tf = ismember(g.labelset, g.labels(vv));
                c = [c find(tf)];
            end
        end

        % which edges contain v
        %  elist = g.edges(v)
        function e = edges(g, v)
        %PGraph.edges Find edges given vertex
        %
        % E = G.edges(V) return the id of all edges from node id V.
            e = [find(g.edgelist(1,:) == v) find(g.edgelist(2,:) == v)];
        end

        function v = vertices(g, e)
        %PGraph.vertices Find vertices given edge
        %
        % V = G.vertices(E) return the id of the nodes that define edge E.
            v = g.edgelist(:,e);
        end


        function [n,c] = neighbours(g, v)
        %PGraph.neighbours Neighbours of a node
        %
        % N = G.neighbours(V) return a vector of ids for all nodes which are
        % directly connected neighbours of node id V.
        %
        % [N,C] = G.neighbours(V) return a vector N of ids for all nodes which are
        % directly connected neighbours of node id V.  The elements of C are the
        % edge costs of the paths to the corresponding node ids in N.
            e = g.edges(v);
            n = g.edgelist(:,e);
            n = n(:)';
            n(n==v) = [];   % remove references to self
            if nargout > 1
                c = g.cost(e);
            end
        end

        function d = cost(g, e)
        %PGraph.cost Cost of edge
        %
        % C = G.cost(E) return cost of edge id E.
            d = g.edgelen(e);
        end

        function p = coord(g, v)
        %PGraph.coord Coordinate of node
        %
        % X = G.coord(V) return coordinate vector, Dx1, of node id V.
            p = g.vertexlist(:,v);
        end

        function showComponent(g, c)
        %PGraph.showcomponent
        %
        % G.showcomponent(C) plots the nodes that belong to graph component C.
            k = g.labels == c;
            showVertices(g, k);
        end


        function c = connectivity(g)
        %PGraph.connectivity Graph connectivity
        %
        % C = G.connectivity() returns the total number of edges in the graph.
            for k=1:g.n
                c(k) = length(g.edges(k));
            end
        end

        function plot(g, varargin)
        %PGraph.plot Plot the graph
        %
        % G.plot(OPT) plot the graph in the current figure.  Nodes
        % are shown as colored circles.
        %
        % Options::
        %  'labels'              Display node id (default false)
        %  'edges'               Display edges (default true)
        %  'edgelabels'          Display edge id (default false)
        %  'MarkerSize',S        Size of node circle
        %  'MarkerFaceColor',C   Node circle color
        %  'MarkerEdgeColor',C   Node circle edge color
        %  'componentcolor'      Node color is a function of graph component

            colorlist = 'bgmyc';

            % show vertices
            holdon = ishold;
            hold on

            % parse options
            opt.componentcolor = false;
            opt.labels = false;
            opt.edges = true;
            opt.edgelabels = false;
            opt.MarkerSize = 8;
            opt.MarkerFaceColor = 'b';
            opt.MarkerEdgeColor = 'b';

            [opt,args] = tb_optparse(opt, varargin);
            
            % set default color if none specified
            if ~isempty(args)
                mcolor = args{1};
            else
                mcolor = 'b';
            end

            % show the nodes as filled circles
            for i=1:g.n
                % for each node
                if opt.componentcolor
                    j = mod( g.component(i)-1, length(colorlist) ) + 1;
                    c = colorlist(j);
                else
                    c = mcolor;
                end
                args = {'LineStyle', 'None', ...
                    'Marker', 'o', ...
                    'MarkerFaceColor', opt.MarkerFaceColor, ...
                    'MarkerSize', opt.MarkerSize, ...
                    'MarkerEdgeColor', opt.MarkerEdgeColor };
                if g.ndims == 3
                    plot3(g.vertexlist(1,i), g.vertexlist(2,i), g.vertexlist(3,i), args{:});
                else
                    plot(g.vertexlist(1,i), g.vertexlist(2,i), args{:});
                end
            end
            % show edges
            if opt.edges
                for e=g.edgelist
                    v1 = g.vertexlist(:,e(1));
                    v2 = g.vertexlist(:,e(2));
                    if g.ndims == 3
                        plot3([v1(1) v2(1)], [v1(2) v2(2)], [v1(3) v2(3)],'k');
                    else
                        plot([v1(1) v2(1)], [v1(2) v2(2)], 'k');
                    end
                end
            end
            % show the edge labels
            if opt.edgelabels
                for i=1:numcols(g.edgelist)
                    e = g.edgelist(:,i);
                    v1 = g.vertexlist(:,e(1));
                    v2 = g.vertexlist(:,e(2));
                    
                    text('String', sprintf('  %g', g.cost(i)), ...
                        'Position', (v1 + v2)/2, ...
                        'HorizontalAlignment', 'left', ...
                        'VerticalAlignment', 'middle', ...
                        'FontUnits', 'pixels', ...
                        'FontSize', 12, ...
                        'Color', 'k');
                end
            end
            % show the labels
            if opt.labels
                for i=1:numcols(g.vertexlist)
                    text('String', sprintf('  %d', i), ...
                        'Position', g.vertexlist(:,i), ... 
                        'HorizontalAlignment', 'left', ... 
                        'VerticalAlignment', 'middle', ... 
                        'FontUnits', 'pixels', ... 
                        'FontSize', 16, ... 
                        'Color', 'b');
                end
            end
            if ~holdon
                hold off
            end
        end

        function v = pick(g)
        %PGraph.pick Graphically select a node
        %
        % V = G.pick() returns the id of the node closest to the point clicked
        % by user on a plot of the graph.
        %
        % See also PGraph.plot.
            [x,y] = ginput(1);
            v = g.closest([x; y]);
        end

        function goal(g, vg)
        %PGraph.goal Set goal node
        %
        % G.goal(VG) for least-cost path through graph set the goal node.  The cost
        % of reaching every node in the graph connected to VG is computed.
        %
        % See also PGraph.path.
            % cost is total distance from goal
            g.goaldist = Inf*ones(1, numcols(g.vertexlist));

            g.goaldist(vg) = 0;
            g.descend(vg);
        end


        function p = path(g, v)
        %PGraph.path Find path to goal node
        %
        % P = G.path(VS) return a vector of node ids that form a path from
        % the starting node VS to the previously specified goal.  The path
        % includes the start and goal node id.
        %
        % See also PGraph.goal.
            p = [v];

            while g.goaldist(v) ~= 0
                v = g.next(v);
                p = [p v];
            end
        end

        function vn = next(g, v)
        %PGraph.next Find next node toward goal
        %
        % V = G.next(VS) return the id of a node connected to node id VS
        % that is closer to the goal.
        %
        % See also PGraph.goal, PGraph.path.
            n = g.neighbours(v);
            [mn,k] = min( g.goaldist(n) );
            vn = n(k);
        end
        
        
        function d = distance(g, v1, v2)
        %PGraph.distance Distance between nodes
        %
        % D = G.distance(V1, V2) return the geometric distance between
        % the nodes with id V1 and V2.
            
            d = g.distance_metric( g.vertexlist(:,v1), g.vertexlist(:,v2));
            
        end

        function [d,k] = distances(g, p)
        %PGraph.distances  Distance to all nodes
        %
        % D = G.distances(V) returns vector of geometric distance from node 
        % id V to every other node (including V) sorted into increasing order
        % by D.
        %
        % [D,W] = G.distances(V) returns vector of geometric distance from node 
        % id V to every other node (including V) sorted into increasing order
        % by D where elements of W are the corresponding node id.
            
            d = g.distance_metric(p, g.vertexlist);
            [d,k] = sort(d, 'ascend');
        end

        function [c,dn] = closest(g, p)
        %PGraph.closest Find closest node
        %
        % V = G.closest(X) return id of node geometrically closest to coordinate X.
        %
        % [V,D] = G.CLOSEST(X) return id of node geometrically closest to coordinate X, and
        % the distance D.
            d = g.distance_metric(p(:), g.vertexlist);
            [mn,c] = min(d);

            if nargin > 1
                dn = mn;
            end
        end

        function display(g)
        %PGraph.display Display state of the graph
        %
        % G.display() displays a compact human readable representation of the
        % state of the graph including the number of vertices, edges and components.
        %
        % See also PGraph.char.
            loose = strcmp( get(0, 'FormatSpacing'), 'loose');
            if loose
                disp(' ');
            end
            disp([inputname(1), ' = '])
            disp( char(g) );
        end % display()

        function s = char(g)
        %PGraph.char Convert graph to string
        %
        % S = G.char() returns a compact human readable representation of the
        % state of the graph including the number of vertices, edges and components.
            s = '';
            s = strvcat(s, sprintf('  %d dimensions', g.ndims));
            s = strvcat(s, sprintf('  %d vertices', g.n));
            s = strvcat(s, sprintf('  %d edges', numcols(g.edgelist)));
            s = strvcat(s, sprintf('  %d components', g.ncomponents));
        end

        %% these are problematic, dont advertise them
        %
        % removing an edge may divide the graph into 2 components, this is expensive
        % to check and currently not implemented
        function delete_edge(g, e)
            g.edgelist(:,e) = [];
            % really need to check if the two halves are connected, is expensive
            % could use path planner
        end

        function delete_node(g, v)
            el = g.edges(v);
            el
            % make the column invalid, really should remove it but this
            % requires changing all the edgelist entries, and the vertex
            % numbers will change...
            g.vertexlist(:,v) = [NaN; NaN];
            g.delete_edge(el);
            g.n = g.n - 1;
        end

        
        function showVertex(g, v)
            %PGraph.showVertex Highlight a vertex
            %
            % G.showVertex(V) highlights the vertex V with a yellow marker.
            
            % TODO allow the line style to be set
            plot(g.vertexlist(1,v), g.vertexlist(2,v), ...
                'LineStyle', 'None', ...
                'Marker', 'o', ...
                'MarkerSize', 12, ...
                'MarkerFaceColor', 'y', ...
                'MarkerEdgeColor', 'y');
        end
        
    end % method

    methods (Access='protected')
        % depth first
        function descend(g, vg)

            % get neighbours and their distance
            for nc = g.neighbours2(vg);
                vn = nc(1);
                d = nc(2);
                newcost = g.goaldist(vg) + d;
                if isinf(g.goaldist(vn))
                    % no cost yet assigned, give it this one
                    g.goaldist(vn) = newcost;
                    %fprintf('1: cost %d <- %f\n', vn, newcost);
                    descend(g, vn);
                else
                    % it already has a cost
                    if g.goaldist(vn) <= newcost
                        continue;
                    else
                        g.goaldist(vn) = newcost;
                        %fprintf('2: cost %d <- %f\n', vn, newcost);
                        descend(g, vn);
                    end
                end
            end
        end

        % breadth first
        function descend2(g, vg)

            % get neighbours and their distance
            for vn = g.neighbours2(vg);
                vn = nc(1);
                d = nc(2);
                newcost = g.goaldist(vg) + d;
                if isinf(g.goaldist(vn))
                    % no cost yet assigned, give it this one
                    g.goaldist(vn) = newcost;
                    fprintf('1: cost %d <- %f\n', vn, newcost);
                    descend(g, vn);
                elseif g.goaldist(vn) > newcost
                    % it already has a cost
                        g.goaldist(vn) = newcost;
                end
            end
            for vn = g.neighbours(vg);
                descend(g, vn);
            end
        end
        function showVertices(g, v)
            for vv=v
                showVertex(g, vv);
            end
        end


        function l = newlabel(g)
            g.curLabel = g.curLabel + 1;
            l = g.curLabel;
            g.ncomponents = g.ncomponents + 1;
            g.labelset = union(g.labelset, l);
        end

        % merge label1 and label2, label2 dominates
        function merge(g, l1, l2)
            ldom = min(l1, l2);
            lsub = max(l1, l2);
            g.labels(g.labels==lsub) = ldom;
            g.ncomponents = g.ncomponents - 1;
            g.labelset = setdiff(g.labelset, lsub);
        end

        function nc = neighbours2(g, v)
            e = g.edges(v);
            n = g.edgelist(:,e);
            n = n(:)';
            n(n==v) = [];   % remove references to self
            c = g.cost(e);
            nc = [n; c];
        end
        
        function d = distance_metric(g, x1, x2)

            % distance between coordinates x1 and x2 using the relevant metric
            % x2 can be multiple points represented by multiple columns
            switch g.measure
                case 'Euclidean'
                   d = colnorm( bsxfun(@minus, x1, x2) );

                case 'SE2'
                    d = bsxfun(@minus, x1, x2);
                    d(3,:) = angdiff(d(3,:));
                    d = colnorm( d );
                otherwise
                    error('unknown distance measure', g.measure);
            end
        end
        
    end
end % classdef
