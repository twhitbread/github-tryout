clear all
Rsun = 696000;
tau = 365.25;
invtau = 1/tau;
eta = 43200000;
maxl = 128;
theta = 0:pi/180:pi;
B0 = 10;
v0 = 1382.4;
Binit = B0*sign(pi/2-theta).*exp(-exp(pi-2.*abs(pi/2-theta)).*v0/(4*eta).*(sin(2.*abs(pi/2-theta))+cos(2.*abs(pi/2-theta))));
alth = zeros(maxl+1,181);
Mll = zeros(maxl);
Yl = zeros(maxl+1,181);

% Store Legendre polynomials and calculate al(0)

leg = zeros(maxl+1,181);
for l=0:maxl
     m = legendre(l,cos(theta));
     leg(l+1,:) = m(1,:);
     Yl(l+1,:) = sqrt((2*l+1)/(4*pi))*leg(l+1,:);
     alth(l+1,:) = Yl(l+1,:).*Binit.*sin(theta);
end

ainit = trapz(theta,alth');

% Velocity profile and derivative

v = -v0.*sin(2*theta).*exp(pi-2.*abs(pi/2 - theta));
deriv = zeros(1,181);
for i = 1:181
if theta(i) <= pi/2
   deriv(i) = -2*v0.*exp(2*theta(i)).*(cos(2*theta(i)) + sin(2*theta(i)));
else
   deriv(i) = 2*v0.*exp(2*(pi-theta(i))).*(sin(2*theta(i))-cos(2*theta(i)));
end
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

odematrix = @(t,a) Mll*a;
tspan = 0:12730;
[T,A] = ode45(odematrix,tspan,ainit);

% Reconstruct B

B = zeros(length(A),181);
B = A*Yl;

h = surf(B);
set(h,'edgecolor','none');
