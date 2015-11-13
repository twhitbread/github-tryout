clearvars -except slnew bip_day
Rsun = 696000;
tau = 365.25;
invtau = 1/tau;
eta = 43200000;
maxl = 128;
theta = 0:pi/180:pi;
B0 = 10;
v0 = 1382.4;
Binit = B0*sign(pi/2-theta).*exp(-exp(pi-2.*abs(pi/2-theta)).*v0/(4*eta).*(sin(2.*abs(pi/2-theta))+cos(2.*abs(pi/2-theta))));

% Store Legendre polynomials and calculate al(0)

leg = zeros(maxl+1,length(theta));
Yl = zeros(maxl+1,length(theta));
alth = zeros(maxl+1,length(theta));

for l=0:maxl
     m = legendre(l,cos(theta));
     leg(l+1,:) = m(1,:);
     Yl(l+1,:) = sqrt((2*l+1)/(4*pi))*leg(l+1,:);
     alth(l+1,:) = Yl(l+1,:).*Binit.*sin(theta);
end

ainit = trapz(theta,alth');

% Velocity profile and derivative

v = -v0.*sin(2*theta).*exp(pi-2.*abs(pi/2 - theta));
deriv = zeros(1,length(theta));
for i = 1:length(theta)
if theta(i) <= pi/2
   deriv(i) = -2*v0.*exp(2*theta(i)).*(cos(2*theta(i)) + sin(2*theta(i)));
else
   deriv(i) = 2*v0.*exp(2*(pi-theta(i))).*(sin(2*theta(i))-cos(2*theta(i)));
end
end

% Block matrices

Mll = zeros(maxl);

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

% Solve for al(t), including new BMRs

tspan = 0:13000;
A2(1,:) = ainit;
sz2 = size(A2);
A = [];

[T1,A1] = ode45(@odematrix,1:219,A2(sz2(1),:),[],Mll); 
sz1 = size(A1);
[T2,A2] = ode45(@odematrix2,219:220,A1(sz1(1),:),[],Mll,slnew(1,:));
sz2 = size(A2);
A = [A;A1;A2];
iprev = 220; 

for i = 219:4240
    for j = 2:5474
        if bip_day(j) == i
            if bip_day(j) == bip_day(j+1)
                slnew(j+1,:) = slnew(j+1,:) + slnew(j,:);
                continue
            elseif bip_day(j)-bip_day(j-1) > 1
                [T1,A1] = ode45(@odematrix,tspan(iprev):tspan(i+1),A2(sz2(1),:),[],Mll); 
                sz1 = size(A1);
                [T2,A2] = ode45(@odematrix2,tspan(i+1):tspan(i+2),A1(sz1(1),:),[],Mll,slnew(j,:));
                sz2 = size(A2);
                A = [A;A1;A2];
                iprev = i+2;
            else
                sza = size(A);
                [T2,A2] = ode45(@odematrix2,tspan(i+1):tspan(i+2),A(sza(1),:),[],Mll,slnew(j,:));
                sz2 = size(A2);
                A = [A;A2];
                iprev = i+2;
            end
        end
    end
end

% Reconstruct B

B = A*Yl;

pcolor(B')
shading interp
colorbar
