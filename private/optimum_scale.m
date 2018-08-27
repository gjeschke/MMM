function [sc,rmsd]=optimum_scale(y,x)
% Optimum scaling of data y_i, so that the root mean square deviation
% rmsd = sum((sc*y_i-x_i)^2) is minimized
% vectors x and y must have the same length
%
% y     vector of data to be scaled
% x     vector of data that should be fitted
% sc    scaling factor
% rmsd  minimized root mean square deviation
%
% G. Jeschke, 2010

sc=sum(x.*y)/sum(y.^2);
diff=sc*y-x;
rmsd=sqrt(sum(diff.^2)/(length(diff)-1));

