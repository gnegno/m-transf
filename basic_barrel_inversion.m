clear all; close all; clc

%% Input
I = checkerboard(50);

[m,n] = size(I);
[X,Y] = meshgrid(1:m,1:n);
xt = X(:) - m/2; 
yt = Y(:) - n/2;

%% From original to barrel

%% Coords handling (1/2)
[theta,r] = cart2pol(xt,yt);

%% Transformation
a = .0001;
s = r + a*r.^3;
[ut1,vt1] = pol2cart(theta, s);

%% Coords handling (2/2)

u1 = reshape(ut1, size(X)) + 200;
v1 = reshape(vt1, size(Y)) + 200;
tmap_B = cat(3,u1,v1);

%resampling
Ib = tformarray(I, [], makeresampler('linear','fill'), [2 1], [1 2], [], tmap_B, .3);


%% From barrel back to original
% 
% %% Coords handling (1/2)
% [theta,r] = cart2pol(xt,yt);
% 
% %% Transformation
% num1 = ((9 * a^2 * r + sqrt(3) * sqrt(27*a^4*r.^2 + 4*a^3)).^(1/3));
% den1 = 2^(1/3) * 3^(2/3) * a;
% 
% num2 = -(2/3)^(1/3);
% den2 = (9 * a^2 * r + sqrt(3) * sqrt(27 * a^4 * r.^2 + 4*a^3)).^(1/3);
% 
% s = num1./den1 + num2./den2;
% [ut1,vt1] = pol2cart(theta, s);
% 
% %% Coords handling (2/2)
% 
% u1 = reshape(ut1, size(X)) + 200;
% v1 = reshape(vt1, size(Y)) + 200;
% tmap_B = cat(3,u1,v1);
% 
% %resampling
% I_ = tformarray(Ib, [], makeresampler('linear','fill'), [2 1], [1 2], [], tmap_B, .3);

%% Visualization
figure;
subplot(131); imshow(I,[]); title('Original');
subplot(132); imshow(Ib,[]); title('Barrel transformation')
% subplot(133); imshow(I_,[]); title('Reconstructed')
