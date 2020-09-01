%       AMM_mainThetaSweep
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
% Usage: AMM_mainThetaSweep
%
% Notes: 
%        Plane waves are: Exp(1i*Kx*r'), Kx = [kx,ky,kz], r = [x,y,z].
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

Period = 1.;                              % [microns, or arbitrary units], all lengths are normalized to Period
lambda = .5/Period;           % Vacuum wavelength, then normalized
theta = 30;                               % Angle of incidence [degrees] measured from the z-axis
TMflag = true;                        % Set problem type (TE, or TM)

nk = TableNK('BK7',lambda*Period);
n1 = 1.0;                                   % superstrate index (real)
n3 = nk(1)+1i*nk(2);             % substrate index (real or complex)

N = 9;                                      % number of Roots to include, integer
FN = 13;                                  % total number of Fourier Modes (odd integer) to expand eigenfunctions by
pathN = 7;                              % for display of roots: number of points/side of 'box' to plot
mpN = 34;                               % Multi-precision decimal digits, when needed

%% ================== Define properties for each grating layer here===============
                    % Grating properties here. Each layer is defined by a
                    % thickness, and a list of index of refractions and
                    % relative widts: [0, x2, x3, . . ., 1]. Note that the
                    % list must start a 0, and end with 1, with xi(n) <= xi(n+1); 
                    % recall that the Period is normalized thus x ranges from 0:1
                    % For each layer length(ei) + 1 = length(xi). 
                    % Note that these conditions are not generally checked in the code below.
% Example of a L layer approximatin to a blaze profile
BlazeAng = 30;		% Blaze angle, degrees
BlazeLay = 11;		  % Number of layers
h = cosd(BlazeAng)*sind(BlazeAng)/Period;			% note that this is in normalized distance (Period == 1)
ei = cell(BlazeLay,1);
xi = cell(BlazeLay,1);
di = cell(BlazeLay,1);
dlayers = zeros(BlazeLay,1);
tn1 = tand(BlazeAng);
tn2 = tand(90 - BlazeAng);
for ii = 1:BlazeLay
	dlayers(ii) = h/BlazeLay;             % vector of grating layer thicknesses
	d1 = cosd(BlazeAng)^2;
	d2 = 1 - d1;
	lh = h*(BlazeLay - ii+.5)/BlazeLay;
	x1 = lh/tn1; %(ii-.5)*h*h/d1/BlazeLay;
	x2 = 1 - lh/tn2; % ((ii-.5)*h*h/d2/BlazeLay);
	ei{ii} = [n1,n3,n1].^2;              % Cell array for list of permittivites (index^2) of each layer (top to bottom)
											  % simple L-layer approximatin to a blaze profile
	xi{ii} = [0., x1, x2,  1.];             % Cell array of list of transition points for each layer
end                        
%% ================== Derived parameters here ================================
FN = FN + mod(FN+1,2);  % enforce odd number of total Fourier basis functions
fN = (FN-1)/2;				   % Fourier-Rayleigh modes (-fN:fN)
k0 = 2*pi/lambda;
kx0 = k0*n1*sind(theta);   
for j = 1:length(xi)
    di{j} = diff(xi{j});         % grating section widths 
end
% This structure holds the grating problem definition
Gparams.k0 = k0;
Gparams.kx0 = kx0;
Gparams.TMflag = TMflag;
Gparams.mpN = mpN;

%% =========================== Theta sweep =================================

% For a theta sweep, the incident k-vector must be evaluated for each update. 

J = complex(zeros(FN,N));     % Fourier-Rayleigh Eigenfunction coefficients for layer material
K = J;                                                 % Fourier-Rayleigh Eigenfunction coefficients for layer material (TM)
Jo =[];
Ko = [];
d2 = [];
mu2 = [];

Roots = complex(zeros(length(dlayers),N));  % array to hold roots for each layer
fnVal = Roots;                        % this array holds the value of the eigenvalue function at the roots (should be close to 0)
itts = Roots;                            % number of function calls needed to find the root
maxAB = -100;					  % for graphing purposes, find maximum of a, b (eigenfunction) coefficients

BC = complex(zeros(FN*2*length(dlayers),N*2*length(dlayers)));      % matrix to hold Bondary Condition Matrix

sweep = linspace(-85,85,590);       %note that theta=90 is not allowed

Bragg = asin(n1*lambda)*180/pi;
DEref = zeros(length(sweep), 5);  % someplace to put the DE, look at 3 modes
DEtrn = DEref;
Transmission = zeros(length(sweep),1);
Reflection = Transmission;
MPneeded = ones(length(sweep),1);		% flags use of MP to find roots

pltIt = false;        % flag to compute and plot eigenvalue equation with roots, and search boxes
kk = 1;                  % easy loop counter
% interest in plotting out the Eigenfunctions in real space (in % x-direction)
xx = linspace(0,1,501)';
EE = complex(zeros(length(xx),N));


% main Theta loop here
for theta =  sweep
     % compute the k-vector components
	while true
		k0 = 2*pi/lambda;
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
		if any(abs(kz1)<1E-5) || any(abs(kz3)<1e-5)
			theta = theta + .001;
		else
			break;
		end
	end

     Gparams.kx0 = kx0;            % keep the grating definition up to date. 
	Gparams.kx0 = kx0;
      % loop through each layer, create/update Stack matrix
      for j = 1:length(dlayers)   % number of layers
  
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
           Roots(j,:) = rts(1:N);
           fnVal(j,:) = fnV(1:N);
           itts(j,:) = itt(1:N);
           
           % update the Fourier-Rayleigh expansion coefficients for the
           % eigenfunctions
           for jj = 1:N                 % for each root find Fourier-Rayleigh eigenfunction expansion
              mu = rts(jj);
		    if ~rtsMPflg(jj), mu = double(mu); end
              [fY,eigCoe] = EigEq(mu,Gparams);
              [an, kan] = AB2FourCoef2(fN,xi{j},eigCoe,kx0);
              J(:,jj) = double(an);
              K(:,jj) = double(kan);  % K for TM mode, with 1/ei in integral
    		    % uncomment to calculate layer eigenfunctions as a function of xx
		    %
		    %mu = double(rts(jj));
              %[fY,eigCoe] = EigEq(mu,Gparams);
		    %[Q,Y,Qp,Yp,E,ex] = QY2x2(xx,xi{j},kx0,eigCoe);
		    %nm = sqrt(trapz(xx,E.*conj(E)./abs(ex)));
              %EE(:,jj) = (E.*exp(-kx0*1i*xx)/nm);
		 end
		%figure(4);
		%plot(xx*Period,[real(EE(:,5)),imag(EE(:,5))]);
		%ylabel('Field Amplitude');
		%title('Eigenfunction N=5');
		%xlabel('position - x');
		
		 % build boundary condition matrix here
           [BC,d2,mu2] = BoundaryMatrix(BC,j,J,Jo,K,Ko,Roots(j,:),kz1,kz3,n1,n3,dlayers(j),Gparams,j==length(dlayers),d2,mu2);
           Jo = J;		% need previous J, K matrices for multi-region layers
           Ko = K;
           if j==1     % for final BC solution need J and d matrices for top and bottom layers
                J1 = J;
                d1 = diag(exp(1i*Roots(j,:)*k0*dlayers(j)));
           end
           if j==length(dlayers)
                Jf = J;
                df = diag(exp(1i*Roots(j,:)*k0*dlayers(j)));
           end
      end
      
      % finish up DE
	ab = zeros(size(BC,1),1);       % ab layer field coefficients
	ab(fN+1) = 2;
	ab = BC\ab;

	  % diagnostic, plot amplitude of field coefficients for each layer as Image
	  % Look at magnitude of field coefficients for downward (an) and upward (bn) eigenfunctions.
	  % Each pair of columns correspond to a grating layer. The expectation is with
	  % 'convergence' the amplitudes of the eigenfunctions should decrease with
	  % increasing N, orders retained. Good visual indication of which
	  % modes contributed the most to the total field solution
     figure(1)           % look at magnitude of field coefficients
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
	if any(ismp(rts))
		MPneeded(kk) = .9;
	end
     kk = kk+1;

     figure(2);
      plot(sweep*Period,[DEref,DEtrn,MPneeded]); grid;
	xlabel('Theta')
	ylabel('Efficiency');
	title('Diffraction Efficiency');
     drawnow;
end
