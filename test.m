[hall,xall] = hist(eld(:),500);
hsort = fliplr(sort(hall));
beginbar = round(length(hall)*0.05);
xrange = 0;

while ~xrange
  gtbarind = find(hall >= hsort(beginbar));
  diffgtbar = diff(gtbarind);
  Nnon1 = find(diffgtbar ~= 1);
  if isempty(Nnon1)
    beginbar = beginbar + 1;
  elseif length(Nnon1) > 1
    beginbar = beginbar - 1;
  else
    xrange = gtbarind(Nnon1)+1:gtbarind(Nnon1+1)-1;
  end
end


t = auto_thresh(eld,[0,5],[]);
xsearch = xall(xall <= t(2));
hsearch = hall(xall <= t(2));