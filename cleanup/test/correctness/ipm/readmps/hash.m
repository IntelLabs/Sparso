%
%  h=hash(key,n)
%
%  Takes a string key and produces a numeric hash value between 1
%  and n.
%
function h=hash(key,n)
h=0;
for i=1:length(key)
  h=h+(key(i)-0)^2*(i*i*i+111);
  h=mod(h,n);
end
h=h+1;