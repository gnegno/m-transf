clear all; close all; clc

A = imread('peppers.png');
A = A(31:330,1:500,:);
figure
%imshow(A)
%title('Original Image','FontSize',14)

%% Ref image
radius = 300;

I = ones(2*radius+1,2*radius+1);

n_circ = 10;
I_circ = I;
for i = 1:n_circ
    I_circ = insertShape(I_circ,'circle', [radius+1,radius+1,(n_circ + 1 - i) * radius/n_circ], 'LineWidth', 2, 'Color', 'black');
end

I_c_line = I_circ;
angles = (0:pi/8:(7/8 * pi))';
points_on_circle = [radius - radius*cos(angles), radius - radius*sin(angles), ...
    2*radius - (radius - radius*cos(angles)), 2*radius - (radius - radius*sin(angles))];

for i = 1:size(points_on_circle,1)
    I_c_line = insertShape(I_c_line,'Line',points_on_circle(i,:),'Color','black');
end

I_c_line(:,1:radius,:) = 1;
A = I_c_line;

%%

t1 = maketform('custom', 2, 2, @conformalForward1, [], []);
t2 = maketform('custom', 2, 2, @conformalForward2, [], []);

uData = [ -1.25   1.25];  % Bounds for REAL(w)
vData = [  0.75  -0.75];  % Bounds for IMAG(w)
xData = [ -2.4    2.4 ];  % Bounds for REAL(z)
yData = [  2.0   -2.0 ];  % Bounds for IMAG(z)

conformal = maketform('custom', 2, 2, [], @myTransf, []);

B = imtransform( A, conformal, 'cubic', ...
    'UData', uData,'VData', vData,...
    'XData', xData,'YData', yData,...
    'Size', [300 360], 'FillValues', 255 );
figure
imshow(B)
title('Transformed Image','FontSize',14)

figure
axIn = conformalSetupInputAxes(axes);
conformalShowInput(axIn, A, uData, vData)
title('Original Image Superposed on Input Plane','FontSize',14)

figure
axOut = conformalSetupOutputAxes(axes);
conformalShowOutput(axOut, B, xData, yData)
title('Transformed Image Superposed on Output Plane','FontSize',14)