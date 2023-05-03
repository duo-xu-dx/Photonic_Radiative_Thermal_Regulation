function psi_v_p = riccati_bessel_psi_p(v, z)
%RICCATI_BESSEL_PSI_P Computes the derivative of the Riccati-Bessel function of the first kind.
%   psi_v_p = RICCATI_BESSEL_PSI_P(v, z) computes the derivative of the 
%   Riccati-Bessel function of the first kind for the order v and argument z.
    psi_v = riccati_bessel_psi(v, z);
    psi_v_minus_1 = riccati_bessel_psi(v - 1, z);
    psi_v_p = psi_v_minus_1 - v * psi_v ./ z;
end