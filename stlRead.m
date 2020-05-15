%STLREAD Reads STL file
%
% [v, f, n, objname] = stlRead(fileName) reads the STL format file (ASCII or
% binary) and returns: 
%
% V (Mx3)   each row is the 3D coordinate of a vertex
% F (Nx3)   each row is a list of vertex indices that defines a triangular face 
% N (Nx3)   each row is a unit-vector defining the face normal
% OBJNAME   is the name of the STL object (NOT the name of the STL file).
%
% Authors::
% - From MATLAB File Exchange by Pau Mico, https://au.mathworks.com/matlabcentral/fileexchange/51200-stltools
% - Copyright (c) 2015, Pau Mico
% - Copyright (c) 2013, Adam H. Aitkenhead
% - Copyright (c) 2011, Francis Esmonde-White
    
%
% stlGetFormat:	identifies the format of the STL file and returns 'binary'
% or 'ascii'. This file is inspired in the 'READ-stl' file written and
% published by Adam H. Aitkenhead
% (http://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation).
% Copyright (c) 2013, Adam H. Aitkenhead. 
%
% stlReadAscii:	reads an STL file written in ascii format. This file is
% inspired in the 'READ-stl' file written and published by Adam H.
% Aitkenhead
% (http://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation).
% Copyright (c) 2013, Adam H. Aitkenhead
%
% stlReadBinary: reads an STL file written in binary format. This file
% is inspired in the 'READ-stl' file written and published by Adam H.
% Aitkenhead
% (http://www.mathworks.com/matlabcentral/fileexchange/27390-mesh-voxelisation).
% Copyright (c) 2013, Adam H. Aitkenhead stlRead:	uses 'stlGetFormat',
%
% 'stlReadAscii' and 'stlReadBinary' to make STL reading independent of the
% format of the file
%
% stlSlimVerts:	finds and removes duplicated vertices. This function is
% written and published by Francis Esmonde-White as PATCHSLIM
% (http://www.mathworks.com/matlabcentral/fileexchange/29986-patch-slim--patchslim-m-).
% Copyright (c) 2011, Francis Esmonde-White.%
%
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the The MathWorks, Inc. nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
%     * Neither the name of the The Christie NHS Foundation Trust nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.


function [v, f, n, name] = stlRead(fileName)
    %STLREAD reads any STL file not depending on its format
    %V are the vertices
    %F are the faces
    %N are the normals
    %NAME is the name of the STL object (NOT the name of the STL file)
    
    format = stlGetFormat(fileName);
    if strcmp(format,'ascii')
        [v,f,n,name] = stlReadAscii(fileName);
    elseif strcmp(format,'binary')
        [v,f,n,name] = stlReadBinary(fileName);
    end
end


function format = stlGetFormat(fileName)
    %STLGETFORMAT identifies the format of the STL file and returns 'binary' or
    %'ascii'
    
    fid = fopen(fileName);
    
    assert(fid > 0, 'Cant find file %s', fileName);
    
    % Check the file size first, since binary files MUST have a size of 84+(50*n)
    fseek(fid,0,1);         % Go to the end of the file
    fidSIZE = ftell(fid);   % Check the size of the file
    if rem(fidSIZE-84,50) > 0
        format = 'ascii';
    else
        % Files with a size of 84+(50*n), might be either ascii or binary...
        % Read first 80 characters of the file.
        % For an ASCII file, the data should begin immediately (give or take a few
        % blank lines or spaces) and the first word must be 'solid'.
        % For a binary file, the first 80 characters contains the header.
        % It is bad practice to begin the header of a binary file with the word
        % 'solid', so it can be used to identify whether the file is ASCII or
        % binary.
        fseek(fid,0,-1); % go to the beginning of the file
        header = strtrim(char(fread(fid,80,'uchar')')); % trim leading and trailing spaces
        isSolid = strcmp(header(1:min(5,length(header))),'solid'); % take first 5 char
        fseek(fid,-80,1); % go to the end of the file minus 80 characters
        tail = char(fread(fid,80,'uchar')');
        isEndSolid = strfind(tail,'endsolid');
        
        % Double check by reading the last 80 characters of the file.
        % For an ASCII file, the data should end (give or take a few
        % blank lines or spaces) with 'endsolid <object_name>'.
        % If the last 80 characters contains the word 'endsolid' then this
        % confirms that the file is indeed ASCII.
        if isSolid & isEndSolid
            format = 'ascii';
        else
            format = 'binary';
        end
    end
    fclose(fid);
end

function [v, f, n, name] = stlReadAscii(fileName)
    %STLREADASCII reads a STL file written in ASCII format
    %V are the vertices
    %F are the faces
    %N are the normals
    %NAME is the name of the STL object (NOT the name of the STL file)
    
    %======================
    % STL ascii file format
    %======================
    % ASCII STL files have the following structure.  Technically each facet
    % could be any 2D shape, but in practice only triangular facets tend to be
    % used.  The present code ONLY works for meshes composed of triangular
    % facets.
    %
    % solid object_name
    % facet normal x y z
    %   outer loop
    %     vertex x y z
    %     vertex x y z
    %     vertex x y z
    %   endloop
    % endfacet
    %
    % <Repeat for all facets...>
    %
    % endsolid object_name
    
    fid = fopen(fileName);
    cellcontent = textscan(fid,'%s','delimiter','\n'); % read all the file and put content in cells
    content = cellcontent{:}(logical(~strcmp(cellcontent{:},''))); % remove all blank lines
    fclose(fid);
    
    % read the STL name
    line1 = char(content(1));
    if (size(line1,2) >= 7)
        name = line1(7:end);
    else
        name = 'Unnamed Object';
    end
    
    % read the vector normals
    normals = char(content(logical(strncmp(content,'facet normal',12))));
    n = str2num(normals(:,13:end));
    
    % read the vertex coordinates (vertices)
    vertices = char(content(logical(strncmp(content,'vertex',6))));
    v = str2num(vertices(:,7:end));
    nvert = length(vertices); % number of vertices
    nfaces = sum(strcmp(content,'endfacet')); % number of faces
    if (nvert == 3*nfaces)
        f = reshape(1:nvert,[3 nfaces])'; % create faces
    end
    
    % slim the file (delete duplicated vertices)
    [v,f] = stlSlimVerts(v,f);
end

function [v, f, n, name] = stlReadBinary(fileName)
    %STLREADBINARY reads a STL file written in BINARY format
    %V are the vertices
    %F are the faces
    %N are the normals
    %NAME is the name of the STL object (NOT the name of the STL file)
    
    %=======================
    % STL binary file format
    %=======================
    % Binary STL files have an 84 byte header followed by 50-byte records, each
    % describing a single facet of the mesh.  Technically each facet could be
    % any 2D shape, but that would screw up the 50-byte-per-facet structure, so
    % in practice only triangular facets are used.  The present code ONLY works
    % for meshes composed of triangular facets.
    %
    % HEADER:
    % 80 bytes:  Header text
    % 4 bytes:   (int) The number of facets in the STL mesh
    %
    % DATA:
    % 4 bytes:  (float) normal x
    % 4 bytes:  (float) normal y
    % 4 bytes:  (float) normal z
    % 4 bytes:  (float) vertex1 x
    % 4 bytes:  (float) vertex1 y
    % 4 bytes:  (float) vertex1 z
    % 4 bytes:  (float) vertex2 x
    % 4 bytes:  (float) vertex2 y
    % 4 bytes:  (float) vertex2 z
    % 4 bytes:  (float) vertex3 x
    % 4 bytes:  (float) vertex3 y
    % 4 bytes:  (float) vertex3 z
    % 2 bytes:  Padding to make the data for each facet 50-bytes in length
    %   ...and repeat for next facet...
    
    fid = fopen(fileName);
    header = fread(fid,80,'int8'); % reading header's 80 bytes
    name = deblank(native2unicode(header,'ascii')');
    if isempty(name)
        name = 'Unnamed Object'; % no object name in binary files!
    end
    nfaces = fread(fid,1,'int32');  % reading the number of facets in the stl file (next 4 byters)
    nvert = 3*nfaces; % number of vertices
    % reserve memory for vectors (increase the processing speed)
    n = zeros(nfaces,3);
    v = zeros(nvert,3);
    f = zeros(nfaces,3);
    for i = 1 : nfaces % read the data for each facet
        tmp = fread(fid,3*4,'float'); % read coordinates
        n(i,:) = tmp(1:3); % x,y,z components of the facet's normal vector
        v(3*i-2,:) = tmp(4:6); % x,y,z coordinates of vertex 1
        v(3*i-1,:) = tmp(7:9); % x,y,z coordinates of vertex 2
        v(3*i,:) = tmp(10:12); % x,y,z coordinates of vertex 3
        f(i,:) = [3*i-2 3*i-1 3*i]; % face
        fread(fid,1,'int16'); % Move to the start of the next facet (2 bytes of padding)
    end
    fclose(fid);
    % slim the file (delete duplicated vertices)
    [v,f] = stlSlimVerts(v,f);
end

function [vnew, fnew]= stlSlimVerts(v, f)
    % PATCHSLIM removes duplicate vertices in surface meshes.
    %
    % This function finds and removes duplicate vertices.
    %
    % USAGE: [v, f]=patchslim(v, f)
    %
    % Where v is the vertex list and f is the face list specifying vertex
    % connectivity.
    %
    % v contains the vertices for all triangles [3*n x 3].
    % f contains the vertex lists defining each triangle face [n x 3].
    %
    % This will reduce the size of typical v matrix by about a factor of 6.
    %
    % For more information see:
    %  http://www.esmonde-white.com/home/diversions/matlab-program-for-loading-stl-files
    %
    % Francis Esmonde-White, May 2010
    
    if ~exist('v','var')
        error('The vertex list (v) must be specified.');
    end
    if ~exist('f','var')
        error('The vertex connectivity of the triangle faces (f) must be specified.');
    end
    
    [vnew, indexm, indexn] =  unique(v, 'rows');
    fnew = indexn(f);
end

