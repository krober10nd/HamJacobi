clearvars; close all; clc;

x=-1:0.10:1; y=x; z=x;

[xg,yg,zg]=meshgrid(x,y,z);

xp = xg(:); yp = yg(:); zp = zg(:); 

figure; scatter3(xp,yp,zp); 

dims = size(xg); c

inod = 3000; 

k(1) = 1;
k(2) = dims(1); 
k(3) = dims(1)*dims(2); 

ndx = inod; 
vi = mod(ndx-1, k(3)) + 1;
vj = (ndx - vi)/k(3) + 1;
ndx = vi;
kpos = vj; 

vi = mod(ndx-1, k(2)) + 1;
vj = (ndx - vi)/k(2) + 1;
ndx = vi;
jpos = vj; 

vi = mod(ndx-1, k(1)) + 1;
vj = (ndx - vi)/k(1) + 1;
ndx = vi;
ipos = vj;

% ---- gather indices using 4 (6 in 3d) edge stencil
npos(1) = inod;
npos(2) = min(jpos,dims(1))*dims(2) + ipos + (kpos-1)*prod(dims(1:2));
npos(3) = max(jpos-2,1)*dims(2) + ipos + (kpos-1)*prod(dims(1:2));
npos(4) = (jpos-1)*dims(2) + min(ipos+1,dims(2)) + (kpos-1)*prod(dims(1:2));
npos(5) = (jpos-1)*dims(2) + max(ipos-1,1) + (kpos-1)*prod(dims(1:2));
npos(6) = (jpos-1)*dims(2) + ipos + min(kpos,dims(3))*prod(dims(1:2));
npos(7) = (jpos-1)*dims(2) + ipos + max(1,kpos-2)*prod(dims(1:2));

hold on; plot3(xg(ipos,jpos,kpos),yg(ipos,jpos,kpos),zg(ipos,jpos,kpos),'ms'); 
hold on; plot3(xg(npos),yg(npos),zg(npos),'rx');
