function [Roots, fnVal, itts, RootIsMP,rts] = FindRootsGRPF(N,Gparams,pltIt)
%       AMM_mainLambdaSweep
%       Analytic Modal Method for TE and TM diffraction grating
%       efficiency calculation of lamellar gratings.
%
%	   Copyright (c) 2020 David Fluckiger
%       Author: David Fluckiger
%       Date: May 2020, revision: 1.0
%
%       This code is provided to encourage the development of AMM method.
%
% Usage: AMM_mainLambdaSweep
%
% This code is derived from Piotr Kowalczyk
% Project homepage: https://github.com/PioKow/GRPF
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
% and associated documentation files (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the Software is furnished to 
% do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
% 
% FindRootsGRPF recovers the roots of the eigenvalue equation for the grating layer defined in Gparams. 
% The following are passed to this routine:
%   N : minimum number of roots to find 
%   Gparams : grating structure, contains definition of grating problem
%   ni : list of index of refractions (though this can also be extracted from Gparams --)
%   lambda : (scaled) wavelength (also can be extracted from Gparams)
%   pltIt : flag to plot eigenvalue/roots found with the contour plot of
%   the eigenvalue equation.

% check to see if uniform layer (single region), if so roots are found analytically
if length(unique(Gparams.ei)) == 1    % have a uniform layer, roots are deterministic
     fN = (N-1)/2;
     Roots = sqrt(Gparams.ei(1) - ((Gparams.kx0-2*pi*(-fN:fN))/Gparams.k0).^2);
     fnVal = zeros(1,N);
     itts = ones(1,N);
     if pltIt                % plot the roots on top of a contour plot of log(abs(f(mu)+1) 
         figure(pltIt+1);
         Uul = max(real(Roots))+.1-.05i;
         Llr =  -.05+1i*max(imag(Roots))+.1i;
         plotEigEq(Uul+.1-.2i,Llr,Gparams)
         hold on;
         plot(1i*real(Roots)+imag(Roots),'+k');
         hold off;
	end
	fnVal = zeros(N,1);
	itts = zeros(N,1);
	RootIsMP = zeros(N,1);
     return;
end

% create sampling points
ni = sqrt(Gparams.ei);
lambda = 2*pi/Gparams.k0;
nrts = 0;
r = lambda/5;					%  r      : initial mesh step
r1 =max(real(ni))*1.5;	
i2 = max(imag(ni))+N*lambda/2;
xb = -0.537;				   %  xb     : real part range begin 
xe = 1.5*r1+3*r;			%  xe     : real part range end 
yb = -0.2;				  %  yb     : imag part range begin 
ye = i2+3*r;				  %  ye     : imag part range end 
bc = [xb,xe,yb,ye,r];

NewNodesCoord = rectdom(xb,xe,yb,ye,r);
Tol = .1e-3; % accuracy 
ItMax=50; % max number of iterations
NodesMax=500000; % max number of nodes
SkinnyTriangle=3; % skinny triangle definition
[nrts,rts,rtsM] = GRPF(NewNodesCoord, Tol, ItMax, NodesMax, SkinnyTriangle, Gparams,bc);
Uu1 = xe+yb*1i;
while nrts<=1.2*N	% might need some more roots, for look at adjacent regions moving along imaginary axis
	imshft = max(4*lambda,lambda*abs(N-nrts)*.6)+3*r;
	yb = .99*ye;
	ye = yb + imshft;
	NewNodesCoord = rectdom(xb,xe,yb,ye,r);
	[nrts1,rts1,rtsMul] = GRPF(NewNodesCoord, Tol, ItMax, NodesMax, SkinnyTriangle, Gparams, bc);
	nrts = nrts+nrts1;
	rts = [rts;rts1];
	rtsM = [rtsM;rtsMul];	% root multiplicity (should == 1)
end
Ll1 = xb+ye*1i;
ul = lambda*(1-1i)*.5;
lr = lambda*(-1+1i)*.5;
Roots = [];
fnVal = [];
itts = [];
RootIsMP = [];
mpN = Gparams.mpN;

for i = 1:length(rts)
	if imag(rts(i))<0, rts(i) = -rts(i); end
	xtry = rts(i);
	if rtsM(i) == 1
		[xn,feVal,itt,~] = rootYasmin(xtry,xtry+ul,xtry+lr,Gparams,true);
		if isnan(xn)
			[xn,feVal,itt,~] = rootYasmin(mp(xtry,mpN),xtry+ul,xtry+lr,Gparams,true);
		end
		if abs(imag(xn)) < 1E-6   % look out for 'conj' on real axis, take only positive one
			xn = abs(real(xn));
		end
		if ~isnan(xn) 
			if ~isempty(Roots)
				if min(min(abs(Roots+conj(xn)')))>1E-5 && min(abs(Roots-xn))>1E-5 
				    Roots = [Roots,xn];         % accumulate roots
				    fnVal = [fnVal,double(feVal)];
				    itts = [itts,itt];
				    RootIsMP = [RootIsMP,ismp(xn)];
				end
			 else
				Roots = [Roots,xn];             % accumulate roots
				fnVal = [fnVal,double(feVal)];
				itts = [itts,itt];
				RootIsMP = [RootIsMP,ismp(xn)];
			end
		end
	else
		xtry = box2Root(xtry+ul,xtry+lr,Gparams,0.1);  
		[xn1,feVal1,itt1,~] = rootYasmin(xtry(1),ul,lr,Gparams,false);  % find both roots, check nan status as with 1 root
		if isnan(xn1)
			[xn1,feVal1,itt1,~] = rootYasmin(mp(xtry(1)),ul,lr,Gparams,false);
		end
		[xn2,feVal2,itt2,~] = rootYasmin(xtry(2),ul,lr,Gparams,false);
		if isnan(xn2)
			[xn2,feVal2,itt2,~] = rootYasmin(mp(xtry(2)),ul,lr,Gparams,false);
		end
		Roots = [Roots,xn1,xn2];
		fnVal = [fnVal,feVal1,feVal2];
		itts = [itts, itt1, itt2];
		RootIsMP = [RootIsMP, 0, 0];
	end
end		
[~,id] = sort(imag(Roots));     % sort roots on imaginary part
Roots = Roots(id);
fnVal = fnVal(id);
itts = itts(id);
RootIsMP = RootIsMP(id);

% if you plot here, will need to set the plot limits dynamically for
% current problem definition!
Uu1 = 4.1-.2i;
Ll1 = 30i - .2;
Uu1 = double(max(real(Roots))*1.1-.2i);
Ll1 = double(max(imag(Roots))*1.1i+min(real(Roots))*1.1-.05);

if pltIt                % plot the roots on top of a contour plot of log(abs(f(mu)+1) 
    figure(pltIt);
    plotEigEq(Uu1+.1-.2i,Ll1,Gparams)
    hold on;
    plot(1i*real(Roots)+imag(Roots),'*r');
     ylabel('\bfreal(\mu)');
    xlabel('\bfimag(\mu)');
    title('\bfRoot (\mu) diagram')
    %axis([-.5,40,-.5,2.5]);
    hold off;
end
