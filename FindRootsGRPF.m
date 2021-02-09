function [Roots, fnVal, itts, RootIsMP] = FindRootsGRPF(N,Gparams,pltIt)
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

% new method to segemeg by rectangles
Tol = lambda/50;		% accuracy to locate root
ItMax=50;					% max number of iterations
NodesMax=500000;   % max number of nodes
SkinnyTriangle=3;	   % skinny triangle definition
nrts = 0;
rts = [];
rtsM = []; 

[~,id] = sort(imag(ni));
sni = ni(id);				% index in increasing imag parg

curni = 1;				% keep track of which index branch is current
xb = -lambda/10;
xe = 0;
ye = -lambda/10;
while nrts < N+5	    % Grab a few more roots than asked for
	di = 1;
	if curni>length(sni), curni=length(sni); end
	if curni <= length(sni)	% not at the last one yet
		di = sum((imag(sni((curni+1):end)) - imag(sni(curni))) <= lambda);
		yb = max(real(sni(curni+(0:di))))+lambda;
		if curni>1 
			xb = xe-4*Tol;		% catch roots close to border
			ye = -lambda;
		end
		xe = xb+imag(sni(curni+di))+2*lambda;
		curni = curni+di+1;
	else			% at or beyond last sni, should have found some root locations
		er = rts(end);
		yb = imag(er) + lambda;
		xb = xe;
		xe = xe + 2*lambda;
	end	
	r = lambda/(4*(di+1));
	bc = [yb,ye,xb,xe,r];
	nodes = rectdom(bc);
	[nrts1,rts1,rtsM1] = GRPF(nodes, Tol, ItMax, NodesMax, SkinnyTriangle, Gparams, bc);
	nrts = nrts+nrts1;
	rts = [rts;rts1];
	rtsM = [rtsM;rtsM1];	% root multiplicity (should == 1)
end
[~,id] = sort(imag(rts));     % sort roots on imaginary part
rts = rts(id);
rtsM = rtsM(id);

% refine any double roots
id = (rtsM>1);
ia = 1:length(rts);
while any(id)
	nx = ia(id);
	Tol = Tol/10;
	for i = 1:length(nx)	% look at each double root and refine
		xn = rts(nx(i));
		rts(nx(i)) = [];
		rtsM(nx(i)) = [];
		nx = nx-1;
		dx = .9*min(abs(rts-xn));
		r = dx/5;
		bc = [real(xn)-dx,real(xn)+dx,imag(xn)-dx,imag(xn)+dx,r];
		nodes = rectdom(bc);
		[~,rts1,rtsM1] = GRPF(nodes, Tol, ItMax, NodesMax, SkinnyTriangle, Gparams, bc);
		rts = [rts;rts1];
		rtsM = [rtsM;rtsM1];
	end
	[~,id] = sort(imag(rts));     % sort roots on imaginary part
	rts = rts(id);
	rtsM = rtsM(id);
	ia = 1:length(rts);
	id =  (rtsM>1);
end

ul = lambda*(1-1i)*.5;
lr = lambda*(-1+1i)*.5;
Roots = [];
fnVal = [];
itts = [];
RootIsMP = [];
mpN = Gparams.mpN;
 
for i = 1:length(rts)
	%if imag(rts(i))<=0, rts(i) = -rts(i); end
	xtry = rts(i);
	%[xn,feVal,itt,~] = rootYasmin(xtry,xtry+ul,xtry+lr,Gparams,false);
	[xn,feVal,itt,~] = rootNewton(xtry,xtry+ul,xtry+lr,Gparams,false);
	%if isnan(xn)
	%	[xn,feVal,itt,~] = rootYasmin(mp(xtry,2*mpN),xtry+ul,xtry+lr,Gparams,false);
	%end
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
