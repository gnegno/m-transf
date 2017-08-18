clear all; close all; clc

%% Input
I = checkerboard(100,5,5);
[m,n] = size(I); 
[X,Y] = meshgrid(1:m,1:n);
xt = X(:) - m/2; 
yt = Y(:) - n/2;

d = m/2; % DISTANCE FROM WHICH THE PATTERN IS LOOKED AT

%% Circular mask
msk = ones(m,n);
msk = 1 - rgb2gray(insertShape(msk, 'filledCircle',[m/2,n/2,m/2],'Color',[0,0,0]));
I = I.*msk;

%% HELMHOLTZ TRANSFORMATION

%% Coords handling (1/2)
[theta, r] = cart2pol(xt,yt);

%% Transformation
tr = @(x)(m * sin ( atan( 2 * x/m ) / 2));
r_t = tr(r) * (m/2)/tr(m/2);
% We multiply by the degree of image magnification that we have at the original
% radius m/2, in order to fit the whole transformed image in the same size
[ut1,vt1] = pol2cart(theta,r_t);

%% Coords handling (2/2)

u1 = reshape(ut1, size(X)) + m/2;
v1 = reshape(vt1, size(Y)) + n/2;
tmap_B = cat(3,u1,v1);

%% Resampling
Ib = tformarray(I, [], makeresampler('linear','fill'), [2 1], [1 2], [], tmap_B, 0);

%% INVERSE TRANSFORMATION

% We already have transformed polar coordinates. We scaled them before, we
% rescale it
r_tt = r;

%% Transformation
k = r_tt / (2*d);
tr_i = @(x)(d * ((2 * x .* sqrt( 1 - x.^2 )) ./ ( sqrt( 1 - 4 * x.^2 .* (1 - x.^2)))));
r_it = tr_i(k);


[ut1,vt1] = pol2cart(theta,r_it);

%% Coords handling (2/2)

u1 = reshape(ut1, size(X)) + m/2;
v1 = reshape(vt1, size(Y)) + n/2;
tmap_B = cat(3,u1,v1);
tmap_B(isinf(tmap_B)) = 0;

%% Resampling
Ic = tformarray(I, [], makeresampler('linear','fill'), [2 1], [1 2], [], tmap_B, 0);

Id = tformarray(Ib, [], makeresampler('linear','fill'), [2 1], [1 2], [], tmap_B, 0);



%% Visualizations
figure(1);
subplot(221); imshow(I,[]); title('Original');
subplot(222); imshow(Ib,[]); title('Barrel transformation')
subplot(223); imshow(Ic,[]); title('Inverse on original')
subplot(224); imshow(Id,[]); title('Inverse on barrel')