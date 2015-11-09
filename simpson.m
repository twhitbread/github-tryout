clear all
Rsun = 696000;
tau = 86400*365.25;
invtau = 1/tau;
eta = 500;
maxl = 128;
theta = 0:pi/360:pi;
B0 = 10;
v0 = 16/1000;
Binit = B0*sign(pi/2-theta).*exp(-exp(pi-2.*abs(pi/2-theta)).*v0/(4*eta).*(sin(2.*abs(pi/2-theta))+cos(2.*abs(pi/2-theta))));
alth = zeros(maxl+1,361);
Mll = zeros(maxl);
Yl = zeros(maxl+1,361);

% Store Legendre polynomials and calculate al(0)

leg = zeros(maxl+1,361);
for l=0:maxl
     m = legendre(l,cos(theta));
     leg(l+1,:) = m(1,:);
     Yl(l+1,:) = sqrt((2*l+1)/(4*pi))*leg(l+1,:);
     alth(l+1,:) = Yl(l+1,:).*Binit.*sin(theta);
end

ainit = trapz(theta,alth');

% Velocity profile and derivative

v0 = 0.016;
v = -v0.*sin(2*theta).*exp(pi-2.*abs(pi/2 - theta));
if theta <= pi/2
   deriv = -2*v0.*exp(2*theta).*(cos(2*theta) + sin(2*theta));
else
   deriv = 2*v0.*exp(2*(pi-theta)).*(sin(2*theta)-cos(2*theta));
end

% Block matrices

for l=0:maxl
   for lp=0
      q = deriv.*sin(theta).*leg(l+1,:).*leg(lp+1,:) + (lp+1).*v.*cos(theta).*leg(l+1,:).*leg(lp+1,:);
      Cll = trapz(theta,q');
      Mll(l+1,lp+1) = -eta*lp*(lp+1)/(Rsun^2) - (2*pi/Rsun)*sqrt((2*l+1)/(4*pi))*sqrt(1/(4*pi))*Cll - invtau;
   end
end

for l=0:maxl
   for lp=1:maxl
      q = deriv.*sin(theta).*leg(l+1,:).*leg(lp+1,:) + (lp+1).*v.*cos(theta).*leg(l+1,:).*leg(lp+1,:) - lp.*v.*leg(lp,:).*leg(l+1,:);
      Cll = trapz(theta,q');
      Mll(l+1,lp+1) = -eta*lp*(lp+1)/(Rsun^2) - (2*pi/Rsun)*sqrt((2*l+1)/(4*pi))*sqrt(1/(4*pi))*Cll - invtau;
   end
end

% Solve for al(t)

tspan = 0:tau;
[T,A] = ode45(@odematrix,tspan,ainit,[],Mll,maxl);

% Reconstruct B

B = zeros(length(A),361);
B = A*Yl;

%h = surf(B);
%set(h,'edgecolor','none');
