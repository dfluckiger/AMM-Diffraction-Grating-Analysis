%       EigEq
%       Evaluate the eigenvalue equation for the layer defined in Gparams
%       at mu
%       
% Usage:
%	[fY,eigCoe] = EigEq(mu,Gparms)
%		mu:  complex number to evaluate the eigenvalue function at
%		Gparms: structure array that defines various parameters specific
%		to a grating layer
%		Gparams.ei (list of region permittivities)
%	     Gparams.k0  (=2*pi/lambda)
%		Gparams.kx0  (=k0*sind(theta))
%		Gparams.TMflag  (=true for TM mode, false for TE mode)
%		Gparams.di  (list of region widths)
%		Gparams.mpN (number of decimal digits for MP calculations)
% Returns
%		fY eigenfunction value at mu
%		eigCoe: convenient structure needed for derivative calculation,
%		holds all intermediate results
%
%	   Copyright (c) 2020 David Fluckiger
%       Author: David Fluckiger
%       Date: May 2020, revision: 1.0
%
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
% and associated documentation files (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the Software is furnished to 
% do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%
function [fY,eigCoe] = EigEq(mu,Gparams)
%comput the n AB (2x2) matrices for a given mu
% comput the Y(mu) function (trace of product of AB's)

di = Gparams.di;
ei = Gparams.ei;
kx0 = Gparams.kx0;
k0 = Gparams.k0;
TMflag = Gparams.TMflag;

n = length(ei);		% number of regions

if ismp(mu)			% if MP preserve MP precision
    AB = mp(zeros(2,2,n),Gparams.mpN);      % the last AB gives Q,Y at x=1;
    acuAB = mp(zeros(2,2,n+1),Gparams.mpN);
else
    AB = zeros(2,2,n);
    acuAB = zeros(2,2,n+1);
end

en = ones(n,1);
if TMflag
    if n>1
        en = [ei(2:end)./ei(1:end-1),1];
    end
end

beta = k0*sqrt(ei-mu.^2);

r = beta.*di;
s = sin(r);         % matlab sinc is sin(pi x)/(pi x), Mathematica is sin(x)/x
c = cos(r);
sb = s./beta;    % make own 'sinc' to avoid exess pi's
it = find(beta==0);
if ~isempty(it)
    sb(it) = di(it);
end
bs = beta.*s;

%compute the AB matrices, but don't accumulate
Y =  [1,0;0,1];  %Y is accumulation of coefficient matrices
acuAB(:,:,1) = Y;
for i = 1:n
    S = [c(i), -bs(i)*en(i); sb(i), c(i)*en(i)];
    Y = Y*S;  %this is the new AB for the segment
    AB(:,:,i) = S; %don't accumulate for derivative calculation
    acuAB(:,:,i+1) = Y;
end

es = 1;
if TMflag
    es = ei(1)/ei(end);
    if isinf(es)
        es = 1;
    end
end

fY =Y(1,1)+es*Y(2,2)-2*cos(kx0);

% load up eigCoe with intermediate results
eigCoe.AB = AB;
eigCoe.acuAB = acuAB;
eigCoe.r = r;
eigCoe.s = s;
eigCoe.c = c;
eigCoe.sb = sb;
eigCoe.bs = bs;
eigCoe.mu = mu;
eigCoe.di = di;
eigCoe.n = n;
eigCoe.beta = beta;
eigCoe.en = en;
eigCoe.TMflag = TMflag;
eigCoe.es = es;
eigCoe.ei = ei;
eigCoe.k0 = k0;
eigCoe.mpN = Gparams.mpN;
