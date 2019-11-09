clearvars; close all; clc; 
%% test cpp stencil call 

elen = 1; 

[xg,yg,zg]=meshgrid(-100:elen:100,-100:elen:100,-100:elen:100); 

xp = xg(:); yp = yg(:); 
zp = zg(:); 

ffun = xg*0 + 1500; 
ffun(zg > 50) = 1750; 

figure; scatter3(xp,yp,zp,50,ffun(:),'filled'); 
%figure; scatter(xp,yp,50,ffun(:),'ffiled'); 

%dims = [size(ffun),1]; 
dims = size(ffun);
dfdx = 100; 
itmax = sqrt(numel(ffun)); 
% fid = fopen('meshsize.txt','w');
% 
% ffun=ffun(:); 
% for i = 1 : length(ffun); 
%     fprintf(fid,'%f\n',ffun(i)); 
% end
% fclose(fid);
% %

smoothed_field = FastHJ( int32(dims), elen, dfdx, int32(itmax), ffun(:));

%figure; scatter(xp,yp,50,smoothed_field,'filled')
figure; scatter3(xp,yp,zp,50,smoothed_field,'filled')