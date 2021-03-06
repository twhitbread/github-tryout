function [bip_day slnew] = bmr
clear all
delta = 180/pi;
Rsun = 696000;
phi = 0:2*pi/180:2*pi;
theta = 0:pi/180:pi;

% Convert to Cartesian coordinates

[th,ph] = meshgrid(theta,phi);

x = sin(th).*cos(ph);
y = sin(th).*sin(ph);
z = cos(th);

% Store Legendre polynomials

maxl = 128;
leg = zeros(maxl+1,181);
Yl = zeros(maxl+1,181);
Ysph = zeros(maxl+1,181);

%progressbar('Getting bmrs...')

for l=0:maxl
  m = legendre(l,cos(theta));
  leg(l+1,:) = m(1,:);
  Yl(l+1,:) = sqrt((2*l+1)/(4*pi))*leg(l+1,:);
  Ysph(l+1,:) = Yl(l+1,:).*sin(theta);
end

% Import data

bmrs = dlmread('bips.txt');
nbips = bmrs(1,1);
bip_day = zeros(nbips,1);
bip_lon = zeros(nbips,1);
bip_lat = zeros(nbips,1);
bip_sep = zeros(nbips,1);
bip_flux = zeros(nbips,1);
bip_tilt = zeros(nbips,1);

% Calculate BMRs in BMR frame

B = zeros(181);
Bnew = zeros(181,5474);
slnew = zeros(nbips,maxl+1);

for i=1:nbips
bip_day(i) = bmrs(i+1,1);
bip_lon(i) = (pi/180)*bmrs(i+1,2);
bip_lat(i) = (pi/180)*bmrs(i+1,3);
bip_sep(i) = bmrs(i+1,4);
bip_flux(i) = bmrs(i+1,5)/((Rsun^2)*(10^10));
bip_tilt(i) = (pi/180)*bmrs(i+1,6);

xp = cos(bip_lat(i)).*x + sin(bip_lat(i)).*z;
yp = sin(bip_tilt(i)).*sin(bip_lat(i)).*x + cos(bip_tilt(i)).*y - sin(bip_tilt(i)).*cos(bip_lat(i)).*z;
zp = -cos(bip_tilt(i)).*sin(bip_lat(i)).*x + sin(bip_tilt(i)).*y + cos(bip_tilt(i)).*cos(bip_lat(i)).*z;
phip = delta.*atan2(yp,xp);
lambdap = delta.*asin(zp);
B = -(bip_flux(i).*(delta^2)/(sqrt(pi).*(bip_sep(i)^3))).*phip.*exp(-(phip.^2 + 2.*(lambdap.^2))/(2.*(bip_sep(i)^2)));

% Integrate in phi

Bnew(:,i) = trapz(phi,B',2);

% Obtain spherical harmonic time-components

for l=0:maxl
sl2 = Ysph(l+1,:).*Bnew(:,i)';
slnew(i,l+1) = trapz(theta,sl2');
end

%progressbar(i/nbips)
end
end