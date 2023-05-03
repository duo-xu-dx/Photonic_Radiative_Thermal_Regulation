function xi_v = riccati_bessel_xi(v, z)
%RICCATI_BESSEL_XI Computes the Riccati-Bessel function of the second kind.
%   xi_v = RICCATI_BESSEL_XI(v, z) computes the Riccati-Bessel function
%   of the second kind for the order v and argument z.
    xi_v = z .* sqrt(pi./(2*z)) .* besselh(v + 1/2, 1, z);
end