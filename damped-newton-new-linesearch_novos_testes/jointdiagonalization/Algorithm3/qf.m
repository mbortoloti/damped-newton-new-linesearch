function [Q] = qf(A)

     [m, n, N] = size(A);
     if m >= n 
         Q = zeros(m, n, N, class(A));
         R = zeros(n, n, N, class(A));
     else
         Q = zeros(m, m, N, class(A));
         R = zeros(m, n, N, class(A));
     end
     
     for k = 1 : N
         
         [q, r] = qr(A(:, :, k), 0);
         
         s = sign(diag(r));
         s(s == 0) = 1;
         
         Q(:, :, k) = bsxfun(@times, q, s.');
         R(:, :, k) = bsxfun(@times, r, conj(s));
         
     end
end