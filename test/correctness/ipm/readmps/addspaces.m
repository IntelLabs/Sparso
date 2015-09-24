%
% r=addspaces(s,n)
%
% Pads string s with spaces to make its length at least n.
%
function r=addspaces(s,n)
r=s;
while (length(r) < n)
  r=[r ' '];
end