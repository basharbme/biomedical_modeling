function [sumRes, resJ] = BallStickSSDC(x, Avox, bvals, qhat)
 
% Extract the parameters
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);
 
% Synthesize the signals
fibdir = [cos(theta)*sin(phi) sin(theta)*sin(phi) cos(phi)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');

%S0 and d values are squared and then square root is applied to make sure S0 and d
%are always positive. absolute value of f is used to constrain the values
%of f between 0 and 1. No constraints are imposed on theta and phi, since
%the trigonometric conversion automatically takes care of these values. 
S = sqrt(S0^2)*(abs(f)*exp(-bvals*sqrt(diff^2).*(fibdotgrad.^2)) + (1-abs(f))*exp(-bvals*diff));
 
% Compute the sum of square differences
sumRes = sum((Avox - S').^2);


