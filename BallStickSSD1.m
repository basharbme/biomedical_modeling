function [sumRes, resJ] = BallStickSSD1(x, Avox, bvals, qhat)
 
% Extract the parameters
S0 = x(1);
diff = x(2);
f = x(3);
theta = x(4);
phi = x(5);
 
% Synthesize the signals
fibdir = [cos(theta)*sin(phi) sin(theta)*sin(phi) cos(phi)];
fibdotgrad = sum(qhat.*repmat(fibdir, [length(qhat) 1])');
 
S = S0*(f*exp(-bvals*diff.*(fibdotgrad.^2)) + (1-f)*exp(-bvals*diff));

figure;
plot(Avox,'o-');hold on;
plot(S,'r.-');
% Compute the sum of square differences
sumRes = sum((Avox - S').^2);


