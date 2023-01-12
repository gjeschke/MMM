function color=color_grade(k,n)
% k-th color grade in a rainbow starting blue and ending red with n grades 

if isnan(k), k=n; end
k=round(k); % allow for real number
if k<1, k=1; end % limit to range
if k>n, k=n; end

if n<=3
    colmap = flipud(hsv(n));
else
    colmap = jet(n);
end
color = colmap(k,:);