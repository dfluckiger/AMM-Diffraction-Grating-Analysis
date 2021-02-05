%       AMM_mainDCSweep
%       Analytic Modal Method for TE and TM diffraction grating
%       efficiency calculation of lamellar gratings.
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
%       This code is provided to encourage the development of AMM method.
%
% Usage: AMM_mainDCSweep
%
% Notes: 
%        Plane waves are: Exp(1i*Kx.r'), Kx = [kx,ky,kz], r = [x,y,z].
%             The harmonic time is Exp(-1i*w t), and is suppressed.
%       The superstrate (z:: -inf:0) is linear, homogeneous, isotropic,
%             infinite half space, described by a real scalar index
%             of refraction: n1. 
%       The substrate (z:: T:+inf) [T = total grating thickness] is linear,
%            homogeneous, isotropic, infinite half-space described 
%            by a real/complex scalar index of refraction: n3.
%       Index of refractions are given in general as n + k*1i, where k>0
%            for absorbing materials.
%       Lamellar grating consisting of N layers, each layer has piecewise
%            constant index of refraction
%       The 'top' of the grating is at z = 0, at the interface of the
%            superstrate and grating. (+z point 'down' through the grating.)
%       The periodic direction of the grating defines the x-coordinate direction,
%            with right handed coordinate assumed. 
%       For TE, TM the plane of incidence is the x-z plane, thus for 
 %           TE, Ey is the only field component (y direction) of the E-field, and
 %           TM Hy is the only field component (y direction) of the H-field.

%% ================= Global Grating Problem Definition ===============================
%clear all;
%close all;

n1 = 1.0;                                   % superstrate index (real)
Period = .5;                              % [microns, or arbitrary units], all lengths are normalized to Period
lambda = .6328/Period;       % Vacuum wavelength, then normalized
theta = 30;                               % Angle of incidence [degrees] measured from the -z-axis
TMflag = true;                        % Set problem type (TE, or TM)	
										% index of refraction of the substrate (table model)
nk = TableNK('AU',lambda*Period);
%n3 = nk(1)+1i*nk(2);              % substrate index (real or complex)
% for the A. Khavasi, et. al., January 15, 2008 / Vol. 33, No. 2 / OPTICS LETTERS, much studied grating set n3 = 10i;
n3 = 10i;

N = 40;                                      % number of eigen-Roots to include in grating layer, integer
FN = 59;                                    % total number of Fourier Modes (odd integer) to expand eigenfunctions in Rayleigh basis

mpN = 34;                               % Multi-precision decimal digits, when needed; default quad-precision (34 decimal digits)


%% ================== Define properties for each grating layer here===============
                    % Grating properties here. Each layer is defined by a
                    % thickness, and a list of index of refractions and
                    % relative widts: [0, x2, x3, . . ., 1]. Note that the
                    % list must start with 0, and end with 1, with xi(n) <= xi(n+1); 
                    % recall that the Period is normalized thus x ranges from 0:1
                    % For each layer length(ei) + 1 = length(xi). 
                    % Note that these conditions are not generally checked in the code below.
% set up a single layer, with the hooks for multiple layers 
dlayers = [.5]/Period;             % vector of grating layer thicknesses
ei = cell(length(dlayers));
xi = cell(length(dlayers));
di = cell(length(dlayers));
% set a two region layer, note that xi will actually be updated in the loop (Duty Cycle sweep)
ei{1} = [n3,n1].^2;         % Cell array for list of permittivites (index^2) of each region, in each layer (top to bottom)
xi{1} = [0,.5,1.];             % Cell array of list of transition points for each layer
                    % example of 2 layer grating (for example)
                    % dlayers = [.5, .782]/Period;
                    % ei{1} = [1.34+.256i, 1.75, 2.1+2.4i].^2;
                    % xi{1} = [0, .25, .5, .68, 1.];
                    % ei{2} = [1.5].^2;  % toss in a uniform layer
                    % xi{2} = [0., 1.];
                    
%% ================== Derived parameters here ================================
FN = FN + mod(FN+1,2);  % enforce odd number of total Fourier basis functions
fN = (FN-1)/2;				   % Fourier-Rayleigh modes (-fN:fN)
% define various k-vector components
k0 = 2*pi/lambda;
kx0 = k0*n1*sind(theta);    
% set the region widths (as diff(xi)) for each layer. xi is a cell array,
% one cell entry for each layer
for j = 1:length(xi)
    di{j} = diff(xi{j});         % grating section widths 
end

% This structure holds the grating problem definition for a particular
% layer, will be updated in the loop below
%Gparams.ei = ei;		%fill in for each layer
%Gparams.di = di;
Gparams.k0 = k0;
Gparams.kx0 = kx0;
Gparams.TMflag = TMflag;
Gparams.mpN = mpN;

%% =========================== Duty Cycle sweep =================================
% allocate various storage locatoins

J = complex(zeros(FN,N));     % Fourier-Rayleigh block Eigenfunction coefficients for layer material
K = J;                                          % Fourier-Rayleigh block Eigenfunction coefficients for layer material (TM)
Jo =[];
Ko = [];
d2 = [];
mu2 = [];

Roots = complex(zeros(length(dlayers),N));			% array to hold roots for each layer
% some diagnostic arrays:
fnVal = Roots;                        % this array holds the value of the eigenvalue function at the roots (should be close to 0)
itts = Roots;                            % number of function calls needed to find the root
maxAB = -100;					  % for graphing purposes, find maximum of a, b (eigenfunction) coefficients

BC = complex(zeros(FN*2*length(dlayers),N*2*length(dlayers)));      % matrix to hold Bondary Condition Matrix

sweep = linspace(0.01,0.99,500);			% Duty Cycle sweep

Bragg = asin(n1*lambda)*180/pi;		% need to trap the Bragg angle with k1z or k3z may be singular
DEref = zeros(length(sweep), 5);			% someplace to put the DE, look at a few modes, 
DEtrn = DEref;
Transmission = zeros(length(sweep),1);	% total transmission and reflection energies
Reflection = Transmission;
MPneeded = ones(length(sweep),1);		% MP flags use of MP to find roots

pltIt = false;        % flag to compute and plot eigenvalue equation with roots, and search boxes
kk = 1;                  % loop counter
% interest in plotting out the Eigenfunctions in real space (in % x-direction)
xx = linspace(0,1,750)';	
EE = complex(zeros(length(xx),N));

% for DC variation, angle dependencies are fixed (calculate once)
% compute the k-vector components, do not allow theta==0 (exactly) or
% theta=Bragg (exatly). Perturb by 17urad.
if abs(theta)<.001
	theta = 0.001;   %need special code for singular angles, so just avoid it
end
if abs(theta-Bragg) < .001
	theta = Bragg+.001;
end
% setup k-vector components for Rayleigh basis
k1 = k0*n1;
k3 = k0*n3;
kx0 = k0*n1*sind(theta);
kx = kx0 - 2*pi*(-fN:fN);   
kz1 = sqrt(k1^2-kx.^2);
kz3 = sqrt(k3^2-kx.^2);
id = find(imag(kz1)<0);       % with choice of phase, attenuated(evanescent) modes must have imag(kz)>0
kz1(id) = -kz1(id);
id = find(imag(kz3)<0);
kz3(id) = -kz3(id);  

for tau =  sweep		% Duty-cycle sweep
	
      % loop through each layer, create/update Stack matrix
      for j = 1:length(dlayers)   % number of layers
		 % for multiple layers, insert code on how to change DC for each layer
			xi{j} = [0., tau, 1.];			% Run DC on layer 1 from 0 to 1
			di{j} = diff(xi{j});
           
		 % for interst in plotting the eigenvalue map update the pltIt flag
		 % pltIt = 0, no plot, pltIt ~= 0, insert which figure number to plot with
		 %if mod(kk,50 )== 0        % choose some interval to examine/plot eigenvalue equation (root locations)
           %     pltIt = 4;
           %else
           %     pltIt = false;
           %end
           
           Gparams.ei = ei{j};  % update current layer definition
           Gparams.di = di{j};
           
		 [rts, fnV, itt, rtsMPflg] = FindRootsGRPF(N,Gparams,pltIt);
		 %[rts, fnV, itt, rtsMPflg] = FindRootsCWA2(N,Gparams,pltIt);
           Roots(j,:) = double(rts(1:N));		%note this converts to double, keep mp for J, K calculation!
           fnVal(j,:) = fnV(1:N);
           itts(j,:) = itt(1:N);
           
           % update the Fourier-Rayleigh expansion coefficients for the
           % eigenfunctions for this layer
           for jj = 1:N                 % for each root find Fourier-Rayleigh eigenfunction expansion
              mu = rts(jj);
		    if ~rtsMPflg(jj), mu = double(mu); end
              [~,eigCoe] = EigEq(mu,Gparams);
              [an, kan] = AB2FourCoef2(fN,xi{j},eigCoe,kx0);  % if root is an MP number, must calculate an/kan in MP arithmetic, then convert to double
              J(:,jj) = double(an);
              K(:,jj) = double(kan);  % K for TM mode, with 1/ei in integral
			    % calculate eigenfunctions as a function of xx here
			    % (uncomment if you want these)
			    % only need double precision for fields
		    mu = double(rts(jj));
              [~,eigCoe] = EigEq(mu,Gparams);
		    [Q,Y,Qp,Yp,E,ex] = QY2x2(xx,xi{j},kx0,eigCoe);
		    nm = sqrt(trapz(xx,E.*conj(E)./abs(ex)));  %eigenfunction normalization 
              EE(:,jj) = (E.*exp(-kx0*1i*xx)/nm);
		 end
		 		 
		 % build boundary condition matrix here
           [BC,d2,mu2] = BoundaryMatrix(BC,j,J,Jo,K,Ko,Roots(j,:),kz1,kz3,n1,n3,dlayers(j),Gparams,j==length(dlayers),d2,mu2);
           Jo = J;		% need previous layer J, K matrices for multi-layer grating
           Ko = K;
           if j==1     % for final BC solution need to save J and d matrices from top and bottom layers
                J1 = J;
                d1 = diag(exp(1i*Roots(j,:)*k0*dlayers(j)));
           end
           if j==length(dlayers)
                Jf = J;
                df = diag(exp(1i*Roots(j,:)*k0*dlayers(j)));
           end
      end
      
     % finish up and compute DE
	ab = zeros(size(BC,1),1);       % ab layer field coefficients
	ab(fN+1) = 2;
	ab = BC\ab;
 
	  % diagnostic, plot amplitude of field coefficients for each layer as Image
	  % Look at magnitude of field coefficients for downward (an) and upward (bn) eigenfunctions.
	  % Each pair of columns correspond to a grating layer. The expectation is with
	  % 'convergence' the amplitudes of the eigenfunctions should decrease with
	  % increasing N, orders retained. Good visual indication of which
	  % modes contributed the most to the total field solution
     figure(1)
     imagesc(10*log10(abs(reshape(ab,N,2*length(dlayers)))),[0,25]);  colorbar
	ylabel('10log_1_0(|c_j|)');
	title('Coe Mag');
	xlabel('Down(odd), Up(even)');
     maxAB = max(maxAB,max(abs(ab(:))));
	
	% reflection, and transmission field amplitudes
     Rm = J1*ab(1:N) + J1*d1*ab(N+(1:N));
     Rm(fN+1) = Rm(fN+1)-1;
     ro = length(ab)-2*N+(1:N);
     Tm = Jf*df*ab(ro) + Jf*ab(ro+N);
	C = 1;
     if TMflag
           C = (n1/n3)^2;  % weighting for TM mode
	end
	
	 % diffraction efficiencies (amplitude of Poynting vector in |z|  direction)
     DE1 = abs(Rm').^2.*real(kz1./(k1*cosd(theta)));
     DE3 =abs(Tm').^2.*real(C*kz3./(k1*cosd(theta)));
     
     Transmission(kk) = sum(DE3);
     Reflection(kk) = sum(DE1);
     DEref(kk,:) = DE1(fN+1+(-2:2));
     DEtrn(kk,:) = DE3(fN+1+(-2:2));
	% indicate if any roots needed multiprecision, just for diagnostic purposes
	if any(ismp(rts))
		MPneeded(kk) = .9;
	end
     kk = kk+1;

     figure(2);
	% take out MPneeded, if not interested
     plot(sweep,[DEref,DEtrn,MPneeded]); grid;
	xlabel('DC')
	ylabel('Efficiency');
	title('Diffraction Efficiency');
     drawnow;
end
% all done
