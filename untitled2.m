clear all; close all; clc


% PARAMETERS (from Balasubramanian et al., 2002)
a = 0.9;
b = 180;
alpha1 = 0.95;
alpha2 = 0.5;
alpha3 = 0.2;

%% Ref image
radius = 300;

I = ones(2*radius+1,2*radius+1);

n_circ = 10;
I_circ = I;
for i = 1:n_circ
    I_circ = insertShape(I_circ,'circle', [radius+1,radius+1,(n_circ + 1 - i) * radius/n_circ], 'LineWidth', 1, 'Color', 'black');
end


I_c_line = I_circ;
angles = (0:pi/8:(7/8 * pi))';
points_on_circle = [radius - radius*cos(angles), radius - radius*sin(angles), ...
    2*radius - (radius - radius*cos(angles)), 2*radius - (radius - radius*sin(angles))];
                
for i = 1:size(points_on_circle,1)
    I_c_line = insertShape(I_c_line,'Line',points_on_circle(i,:),'Color','black');
end

I = rgb2gray(I_c_line(:,(radius-1):end,:));

%%
x_center = 0;
y_center = radius;

[m,n] = size(I);
[X,Y] = meshgrid(1:n,1:m);
xt1 = X(:) - x_center; 
yt1 = Y(:) - y_center;

[theta,r] = cart2pol(xt1,yt1);

z = r .* exp(1i * alpha1 * theta);
% z_t = 60*log((z+a));

% theta = unwrap(angle(z_t));
% r = abs(z_t);
% 
% [xt2,yt2] = pol2cart(theta,r);
xt2 = real(z_t);
yt2 = imag(z_t);


%%
U = reshape(xt2, size(X)) + x_center;
V = reshape(yt2, size(Y)) + y_center;

Ib = interp2(X,Y,I,U,V);

% tmap_B = cat(3,u1,v1);
% 
% %resampling
% Ib = tformarray(I, [],makeresampler({'nearest','cubic'},'fill'), [2 1], [1 2], [], tmap_B, 0);

figure(2); imshow(I)
figure(3); imshow(Ib)

