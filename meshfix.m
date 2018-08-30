function [faces, vertices] = meshfix(faces, vertices, joinMultipleComponents)
%MESHFIX Summary of this function goes here
%   Detailed explanation goes here

narginchk(2,3);
nargoutchk(0,2);
if nargin < 3 || isempty(joinMultipleComponents)
    joinMultipleComponents = false;
end

validateattributes(faces, {'numeric'}, {'2d', 'ncols', 3, 'integer', 'positive'}, ...
    mfilename, 'faces', 1);
validateattributes(vertices, {'numeric'}, {'2d', 'ncols', 3, 'finite', 'real'}, ...
    mfilename, 'vertices', 2);
validateattributes(joinMultipleComponents, {'numeric', 'logical'}, ...
    {'scalar', 'binary'}, mfilename, 'joinMultipleComponents', 3);

if exist('mat_meshfix', 'file') ~= 3
    fprintf('Compiling mat_meshfix\n');
    prevDir = pwd;
    cd(fileparts(mfilename('fullpath')));
    if isfolder('build')
        rmdir('build', 's');
    end
    status = system('./build_meshfix.sh', '-echo');
    assert(status==0, [mfilename, ':compilationFailed'], ...
        'Compilation of ''mat_meshfix'' failed, check logs');
    cd(prevDir);
end
[faces, vertices] = mat_meshfix(double(faces)-1, double(vertices), joinMultipleComponents);
faces = faces + 1;


end

