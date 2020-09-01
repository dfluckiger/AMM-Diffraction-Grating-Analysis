%       rootYasmin
%       Find a root of the eigvalue equation given the starting location mu
%       
% Usage:
%	[x1,fY,itt,eigCoe] = rootYasmin(mu,Uul,Llr,Gparms,mpN,perFlag)
%		mu:  complex number to evaluate root of eigenvalue function
%		Uul, Llr, are complex numbers that define a search volume (prevent walk-off)
%		Gparms: structure array that defines various parameters specific
%		to a grating layer
%		Gparams.ei (list of region permittivities)
%	     Gparams.k0  (=2*pi/lambda)
%		Gparams.kx0  (=k0*sind(theta))
%		Gparams.TMflag  (=true for TM mode, false for TE mode)
%		Gparams.di  (list of region widths)
%		Gparams.mpN  number of digits to use if MP is used
%		perFlag  flag to employ bound constraints (do not allow root to
%		walk out of box: true, false for turn off check)
% Returns
%		x1 root, of 1E10 if no convergnece
%		fY, value of eigenfunction function at x1 (should be close to zero)
%		itt, number of itterations
%		eigCoe, intermediate results from evaluation of eigenfunction at	x1
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
function [x1,fY,itt,eigCoe] = rootYasmin(mu,Uul,Llr,Gparams,perFlag)

% compute the root nearest mu by Yasmin (2016), 
% 2-point 4th order with derivative, optimal Hermite interpolation
% must remain in box=[Uul,Llr]; real(Uul)>real(LLr), imag(Llr)>imag(Uul),
% enforce periodic boundary conditions so root does not walk off if perFlag==true

delRe = real(Uul-Llr);          %dimensions of containing box
delIm = imag(Llr-Uul);
zRef = real(Llr)+imag(Uul)*1i;  %lower left corner
mpN = Gparams.mpN;

if ismp(mu)         % stopping criteria here, set a 'zero' location when abs(f(z)) is small
    lim2 = 1E-6;
    itlim = round(1.25*mpN);      % limit on number of search steps/itterations
else
    lim2 = 1E-8;
    itlim = 24;         % for 'good' starting value often, only 1 or a few iterations are needed
end                         % however, for poor guess (say center of large box) allow for sufficient passes 
                                % Note: approach to zero is not uniform, and depends highly on local stucture of
                                % f(z), which has lots of wiggles in it (exponential*trigonometric polynomial)
% this enforces boundary conditions, keeping estimate in the
% defined box limits:
xn = mu;   
if perFlag
	xn = zShift(xn,Uul,Llr);
end
% initial calculation of eigenfunction and derivative at current estimate (xn)
[fY,eigCoe] = EigEq(xn,Gparams);
dY = dEigEq(eigCoe);
% initialize function calls, count derivative call same as function call
% (even though it is n times more computationally expensive, n = number of regions)
itt = 2; 
% check to see if root is already found
if abs(fY) < lim2
    x1 = xn;
    return
end

while true  % main itteration loop 
    % find the new root estimate. First step is a Newton-Rasphon step
     stp = fY/dY;     
     yn = xn - stp;
     % apply BC
	if perFlag
		yn = zShift(yn,Uul,Llr);
	end
     % find function value for new root, update count
     [fZ,~] = EigEq(yn,Gparams);     
     itt = itt+1;
     % check to see if EigEq had a problem (this is highly unlikely)
     if isnan(fZ)
         x1 = nan; fY = 1E15; itt = 0;  % if so flag with nan
		% if x1 returns as nan, then invoke MP and try again
		return;
     end
     % check to see if root found, and also if stp is approaching eps,
     % accept root if function is small too.
     if abs(stp)<1E-15 && abs(fZ)<1E-5 || abs(fZ)< lim2
		x1 = yn; 
		fY = fZ;
		return;
	end
     
	sts = fZ/(2*((fY-fZ)/(xn-yn))-dY);  %Yasmin eqn 2.7, based on Hermite Interpolation
     if abs(xn-yn)<8E-16 && abs(fZ-conj(fY))<1E-6 && abs(fZ)<1E-6 	% solution is oscillating around root
		x1 = yn;  fY = fZ;
		return;
	end	    
     
     xn = yn - sts;     %new root estimate, check periodic BC
	if perFlag
		xn = zShift(xn,Uul,Llr);
	end
     % evaluate eigenfunction and derivative here
     [fY,eigCoe] = EigEq(xn,Gparams); 
     itt = itt+1;
     if isnan(fY)
         x1 = nan; fY = 1E10; 
       return;
     end
    if abs(sts)<1E-15 && abs(fY)<1E-5 ||  abs(fY)< lim2
        x1 = xn; 
       return
    end    
    dY = dEigEq(eigCoe);

    % check to see if run out of itterations:
     if ~ismp(mu) && itt>itlim 
         x1 = nan; fY = 1E10;
        return; 
     else
         if ismp(mu) && itt>itlim
             x1 = nan; fY = 1E10;
             return;
         end
     end
        
end
     