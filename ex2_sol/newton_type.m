function W = newton_type(w0, F, M, numit)
% implementation of a Newton type algorithm for finding the root of
% function F(w), i.e., finding w such that F(w) = 0

% Inputs:   w0:     initial guess (column vector)
%           F:      function handle of F
%           M:      function handle of a approximation of the Jacobian of F
%           numit:  number of iterations performed

% Returns:
%           W:      iteration history of w
W = w0;

for k = 1:numit
    w = W(:,k);
    W(:,k+1) = w - M(w) \ F(w);     % Newton-type step
end
