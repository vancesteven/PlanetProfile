function out = Helgerud_sI(P,T)
% input (P,T) to get structure contains Vp, Vs, shear, K and rho

% From Ning et al 2014
%alpha=1/v *dv/dt

out.alpha=((2*3.5038*10^-4)*T+0.2517)./(3.5038*10^-4*(T.^2)+0.2517*T+1610.9373);


% from Helgerud et al. (2009)
T=T-273.15; % temperature in celcius 
out.Vp = -1.84.*T + 0.31.*P + 3766;
out.Vs = -0.892.*T -0.1.*P + 1957;
out.shear = -4.2e-3.*T + 9e-5.*P + 3.541;
out.K = -1.9e-2.*T + 3.8e-3.*P + 8.39;
out.rho = (-2.3815e-4.*T + 1.1843e-4.*P +0.91801)*1e3;

