function mps2mtx(filename);

problem=readmps([filename '.mps']);
mmwrite([filename '-l.mtx'], full(problem.lbnds));
mmwrite([filename '-u.mtx'], full(problem.ubnds));

nslacks = 0;
m = size(problem.A, 1);

pidx = -1;
p = [];
if size(problem.ranges,1) == 0
  for i = 1:size(problem.rowtypes, 2)
    type = problem.rowtypes(i);
    if strcmp(type, 'L')
      nslacks = nslacks + 1;
      problem.A = [ problem.A [ zeros(i - 1, 1); 1; zeros(m - i, 1) ] ];
    elseif strcmp(type, 'G')
      problem.A(i,:) = -problem.A(i,:);
      problem.A = [ problem.A [ zeros(i - 1, 1); 1; zeros(m - i, 1) ] ];
      nslacks = nslacks + 1;
    elseif strcmp(type, 'N')
      if pidx > 0
        'error duplicated objective function'
      else
        p = problem.A(i,:);
        pidx = i;
      end
    end
  end
else
  'error: cannot handle ranges'
end

if pidx > 0
  problem.A(pidx,:) = [];
  problem.rhs(pidx,:) = [];
  p = [p zeros(1, size(problem.A, 2) - size(p, 2))];
end

if size(problem.ubnds,2) > 0
  'error: cannot handle ranges!'
  if size(problem.rhs) == 0
    mmwrite([filename '-b.mtx'], -problem.A*full(problem.ubnds)');
  else
    mmwrite([filename '-b.mtx'], full(problem.rhs) - problem.A*full(problem.ubnds)');
  end
else
  mmwrite([filename '-b.mtx'], full(problem.rhs));
end
mmwrite([filename '-r.mtx'], full(problem.ranges));
mmwrite([filename '-A.mtx'], problem.A);
if pidx > 0
  mmwrite([filename '-p.mtx'], full(p)');
end
