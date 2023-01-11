function ang = subspace_angles(A,B)
%     References
%     ----------
%     .. [1] Knyazev A, Argentati M (2002) Principal Angles between Subspaces
%            in an A-Based Scalar Product: Algorithms and Perturbation
%            Estimates. SIAM J. Sci. Comput. 23:2008-2040.


% Compute orthonormal bases
QA = orth(A);
QB = orth(B);

% svd for cosine
QA_H_QB = QA'* QB;
sigma = svd(QA_H_QB);

% compute matrix B
if size(QA,2) >= size(QB,2)
    B = QB - (QA * QA_H_QB);
else
    B = QA - (QB * QA_H_QB');
end

% svd for sine
mask = sigma.^2  >= 0.5;
if sum(mask) == 0
    ang = zeros(size(mask));
    return
end

if any(mask)
    sigma_B = svd(B);
    sigma_B(sigma_B<-1) = -1;
    sigma_B(sigma_B>1)  = 1;
    mu_arcsin = asin(sigma_B);
else
    mu_arcsin = 0.0;
end

% Compute the principal angles - smallest sigma belongest to largest angle
s = flip(sigma);
s(s<-1) = -1;
s(s>1) = 1;
mu_arccos = acos(s);

% angles are elements of mu_arcsin where mask is true, elements of
% mu_arccos otherwise
ang = nan(size(mask));
ang(mask) = mu_arcsin(mask);
ang(~mask) = mu_arcsin(~mask);


ang = rad2deg(ang);

end




















