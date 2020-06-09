function W = newton_type(w0, F, M, numit)

W = w0;
for k = 1:numit
    w_k = W(:, k);
    
    %add TOL for gradient value
    W(:, k+1) = w_k - M(w_k) \ F(w_k);
end
