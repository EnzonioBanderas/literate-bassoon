function ind_blink = identifyBlinks(d)

  %t = ed.time;
  %d = ed.pupilsize_left;
  md = mean(d);
  %sf = round(1/mean(diff(t)),3);
  d2 = ter_smoothCurve(d,10,'nocomment');
  d3 = [0;diff(d2)];
  ind_ons = find(diff(d==0)== 1);
  ind_end = find(diff(d==0)==-1);
  if d(1) == 0 
    ind_ons = [1;ind_ons];
  end
  if d(end) ==0
    ind_end = [ind_end;numel(d)];
  end

  ind_ons2 = ind_ons;
  for i=1:numel(ind_ons)
    i1 = ind_ons(i);
    while abs(d3(i1))>0.005 *md && i1>1
      i1 = i1-1;
    end
    ind_ons2(i)=i1;
  end
 
  ind_end2 = ind_end;
  for i=1:numel(ind_end)
    i1 = ind_end(i);
    while abs(d3(i1))>0.005 *md && i1<=numel(d)
      i1 = i1+1;
    end
    ind_end2(i)=i1;
  end
  
  ind_blink0  = false(size(d));
  ind_blink2 = false(size(d));
  for i=1:numel(ind_ons)
    ind_blink0(ind_ons(i):ind_end(i))=true;
    ind_blink2(ind_ons2(i):ind_end2(i))=true;
  end
  ind_blink=ind_blink2;
  
%   d4 = d;
%   d5 = d;
%   d4(ind_blink0)=nan;
%   d5(ind_blink2)=nan;
%   clf
%   plot(t,d);
%   hold on;
%   plot(t,d4);
%   plot(t,d5);
%   legend('main','eb1','eb2')
  
end