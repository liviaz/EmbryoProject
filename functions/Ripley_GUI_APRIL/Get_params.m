function [sidex sidey sidez lambda] = Get_params(dists)
% Get_params
% Input: dists - a cell array with distributions
% Output: [sidex sidey lambda] The mean values from dists
%
% License: RipleyGUI is distributed free under the conditions that
% (1) it shall not be incorporated in software that is subsequently sold; 
% (2) the authorship of the software shall be acknowledged in any publication that uses results generated by the software; 
% (3) this notice shall remain in place in each source file. 

sidex = 0;
sidey = 0;
sidez = 0;
lambda = 0;
for i=1:length(dists)
    X = dists{i}.X;
    Y = dists{i}.Y;
    Z = dists{i}.Z;
    xx = max(X) - min(X);
    yy = max(Y) - min(Y);
    zz = max(Z) - min(Z);
    lambdalambda = length(X)/(xx*yy*zz);
    sidex = sidex + xx;
    sidey = sidey + yy;
    sidez = sidez + zz;
    lambda = lambda + lambdalambda;
end
sidex = sidex/length(dists);
sidey = sidey/length(dists);
sidez = sidez/length(dists);
lambda = lambda/length(dists);

end %end Get_params