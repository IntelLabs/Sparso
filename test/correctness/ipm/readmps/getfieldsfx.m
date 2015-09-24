%
% [line,fields]=getfieldsfx(fid)
%
% Reads the next nonempty, non-comment line of the input file and
% returns the fields of this line in a cell array.
%
function [line,fields]=getfieldsfx(fid)
%
% Read a line and parse it.
%
line=fgetl(fid);
if (line == -1)
  error('Unexpected end of file.');
end
line=stripcomments(line);
fields=parsefieldsfx(line);
%
% Skip through blank and comment lines until we get something real.
%
while ((length(line)==0) | (length(fields)==0))
  line=fgetl(fid);
  if (line == -1)
    error('Unexpected end of file.');
  end
  line=stripcomments(line);
  fields=parsefieldsfx(line);
end

