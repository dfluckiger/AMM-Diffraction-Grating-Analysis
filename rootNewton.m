function [xn,fY,itt,eigCoe] = rootNewton(mu,Uul,Llr,Gparams,perFlag)

% compute the root nearest mu by Yasmin (2016) method, with derivative
% 2-point 4th order with derivative, optimal Hermite interpolation
% must remain in box=[Uul,Llr]; real(Uul)>real(LLr), imag(Llr)>imag(Uul),
% enforce periodic boundary conditions so root does not walk off

% auto-mp

% this enforces periodic boundary conditions, keeping estimate in the
% defined box limits:
xn = mu;   
if perFlag
	xn = zShift(xn,Uul,Llr);
end
% initial calculation of eigenfunction and derivative at current estimate (xn)
[fY,eigCoe] = EigEq(xn,Gparams);
dY = dEigEq(eigCoe);
m = ceil(log10(dY));
if m>6
	mpN = 14+m;
	if mpN<34, mpN = 34; end		%take advanatage of quad-double code
	xn = mp(xn,mpN);
end
lim2 = 1E-6;
itlim = 25;      % limit on number of search steps/itterations
% initialize function calls, count derivative call same as function call
itt = 2; 
% check to see if root is already found
if abs(fY) < lim2
    return
end

while true  % main itteration loop 
    % find the new root estimate. First step is a Newton-Rasphon step
     stp = fY/dY;     
     xn = xn - stp;
     % apply periodic BC
	if perFlag
		xn = zShift(xn,Uul,Llr);
	end
     % find function value for new root, update count
     [fY,~] = EigEq(xn,Gparams);     
     itt = itt+1;
     % check to see if EigEq had a problem (this is probably highly unlikely
     if isnan(fY)
         xn = nan; fY = 1E15; itt = 0;  % if so flag with nan
		return;
     end
     % check to see if root found, and also if stp is approaching eps,
     % accept root if function is small too.
	if abs(fY)< lim2 
		return;
	end
     [fY,eigCoe] = EigEq(xn,Gparams); 
     dY = dEigEq(eigCoe);

	% check to see if run out of itterations:
     
    if itt>itlim
	   xn = nan; fY = 1E10;
	   return;
    end
         
end
     