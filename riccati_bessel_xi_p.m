function xi_v_p = riccati_bessel_xi_p(v, z)
%RICCATI_BESSEL_XI_P Computes the derivative of the Riccati-Bessel function of the second kind.
%   xi_v_p = RICCATI_BESSEL_XI_P(v, z) computes the derivative of the 
%   Riccati-Bessel function of the second kind for the order v and argument z.
    xi_v = riccati_bessel_xi(v, z);
    xi_v_minus_1 = riccati_bessel_xi(v - 1, z);
    xi_v_p = xi_v_minus_1 - v * xi_v ./ z;
end