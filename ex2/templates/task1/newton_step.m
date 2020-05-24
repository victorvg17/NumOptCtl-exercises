function w_next = newton_step(hess, grad, w)
    H = hess(w(1), w(2));
    w_next = w - inv(H) * grad(w(1), w(2));
end