function psi_v = riccati_bessel_psi(v, z)
%RICCATI_BESSEL_PSI Computes the Riccati-Bessel function of the first kind.
%   psi_v = RICCATI_BESSEL_PSI(v, z) computes the Riccati-Bessel function
%   of the first kind for the order v and argument z.
    psi_v = z .* sqrt(pi./(2*z)) .* besselj(v + 1/2, z);
end