function U = myTransf( X, ~ )
%MYTRANSF Summary of this function goes here
%   Detailed explanation goes here
Z = complex(X(:,1),X(:,2));
W = log(Z);
U(:,2) = imag(W);
U(:,1) = real(W);

end

