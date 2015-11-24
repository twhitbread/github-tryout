clear all
%clearvars -except slnew bip_day
Rsun = 696000;
tau = 10*24*365;
eta = 1800000;
maxl = 128;
theta = 0:pi/180:pi;
B0 = 10;
v0 = 57.6;
v = -sin(2*theta).*exp(pi-2.*abs(pi/2 - theta));
v0 = v0/max(abs(v));
v = v*v0;
Binit = B0*sign(pi/2-theta).*exp(-exp(pi-2.*abs(pi/2-theta)).*v0/(4*eta).*(sin(2.*abs(pi/2-theta))+cos(2.*abs(pi/2-theta))));

[bip_day slnew] = bmr;

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

deriv = zeros(1,length(theta));
for i = 1:length(theta)
if theta(i) <= pi/2
   deriv(i) = -2*v0.*exp(2*theta(i)).*(cos(2*theta(i)) + sin(2*theta(i)));
else
   deriv(i) = 2*v0.*exp(2*(pi-theta(i))).*(sin(2*theta(i))-cos(2*theta(i)));
end
end

% Block matrices

Cll = zeros(maxl+1);

for l=0:maxl
   for lp=0
      q = (deriv.*sin(theta) + (lp+1).*v.*cos(theta)).*leg(l+1,:).*leg(lp+1,:);
      Cll(l+1,lp+1) = trapz(theta,q');
      Cll(l+1,lp+1) = Cll(l+1,lp+1)*(2*pi/Rsun)*sqrt((2*l+1)/(4*pi))*sqrt((2*lp+1)/(4*pi));
   end
   for lp=1:maxl
      q = (deriv.*sin(theta) + (lp+1).*v.*cos(theta)).*leg(l+1,:).*leg(lp+1,:) - lp.*v.*leg(lp,:).*leg(l+1,:);
      Cll(l+1,lp+1) = trapz(theta,q');
      Cll(l+1,lp+1) = Cll(l+1,lp+1)*(2*pi/Rsun)*sqrt((2*l+1)/(4*pi))*sqrt((2*lp+1)/(4*pi));
   end
   %Mll = diag(-eta*l.*(l+1)/(Rsun^2) - invtau);% - Cll;
end

l=0:maxl;
if (tau > 0)
  Mll=diag((-eta/(Rsun^2))*l.*(l+1) - 1/tau);
else
  Mll=diag((-eta/(Rsun^2))*l.*(l+1));
end
Mll = Mll - Cll;

% Solve for al(t), including new BMRs

% First BMR

A = [];
T = [];

[T1,A1] = ode45(@odematrix,0:5232,ainit,[],Mll);
A = [A;A1];
T = [T,T1];

sz1 = size(A1);
abmr = A1(sz1(1),:)+slnew(1,:);
[T1,A1] = ode45(@odematrix,5232:5544,abmr,[],Mll);
A = [A;A1];
T = [T;T1];
iprev = 231*24;
sza = size(A);

for i = 231:12564
    for j = 2:5473
        if bip_day(j) == i
            if bip_day(j) == bip_day(j+1)
                slnew(j+1,:) = slnew(j+1,:) + slnew(j,:);
                bip_day(j) = bip_day(j-1);
                continue
            elseif bip_day(j)-bip_day(j-1) ~= 1
                [T1,A1] = ode45(@odematrix,iprev:bip_day(j+1)*24,A(sza(1),:)+slnew(j,:),[],Mll); 
                A = [A;A1];
                T = [T;T1];
                sza = size(A);
                iprev = bip_day(j+1)*24; 
            else
                [T1,A1] = ode45(@odematrix,bip_day(j)*24:bip_day(j+1)*24,A(sza(1),:)+slnew(j,:),[],Mll);
                A = [A;A1];
                T = [T;T1];
                sza = size(A);
                iprev = bip_day(j+1)*24;
            end
        end
    end
end

% Final BMR

if iprev == 12730*24;
[T1,A1] = ode45(@odematrix,iprev:13000*24,A(sza(1),:)+slnew(5474,:),[],Mll); 
A = [A;A1];
T = [T;T1];
end

% Reconstruct B

B = fliplr(A*Yl);

pcolor(B')
title('Longitude-Averaged B_r','FontSize',14)
shading interp
colorbar
caxis([-10,10])
xlabel('Time (years)','FontSize',14)
set(gca,'XTick',1:48833:292999)
set(gca, 'XTickLabel',0:5:30)
ylabel('Latitude (degrees)','FontSize',14)
set(gca, 'YTick',[1;91;181])
set(gca, 'YTickLabel', [-90;0;90])

