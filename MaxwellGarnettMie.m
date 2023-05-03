function epsilon_eff = MaxwellGarnettMie(lambda,num_filler,f,eps_host,eps_filler,R)
%MAXWELLGARNETTMIE Computes the effective dielectric constant using a size-dependent extention of the MG formula form the Mie theory.
%   epsilon_eff = MAXWELLGARNETTMIE(lambda, num_filler, f, eps_host, eps_filler, R) calculates the effective 
%   dielectric constant (epsilon_eff) for a composite material using the Maxwell-Garnett effective medium 
%   approximation and Mie scattering. The function takes the following inputs:
%   
%   lambda:     Wavelengths (in the same units as R) as a vector.
%   num_filler: Number of filler materials in the composite.
%   f:          Volume fractions of the filler materials as a vector of length num_filler.
%   eps_host:   Dielectric constants of the host material at each wavelength as a vector.
%   eps_filler: Dielectric constants of the filler materials at each wavelength as a 
%               num_filler x numel(lambda) matrix.
%   R:          Radii of the filler materials as a vector of length num_filler (in the same units as lambda).
%
%   epsilon_eff: The effective dielectric constant of the composite material as a vector with the same length as lambda.
%
%   The function calculates the effective dielectric constant using the Maxwell-Garnett effective medium
%   approximation, which models the composite material as a host material with embedded spherical
%   inclusions (filler materials). It takes into account Mie scattering by using Riccati-Bessel functions
%   to describe the scattering behavior of the spherical inclusions.

RHS_i = zeros(num_filler,numel(lambda));
alpha_1i = zeros(num_filler,numel(lambda));
for i = 1:num_filler
    for j = 1:numel(lambda)
        alpha_1i(i,j) = (sqrt(eps_filler(i,j))*riccati_bessel_psi(1,2*pi*R(i)*sqrt(eps_filler(i,j))/lambda(j))*riccati_bessel_psi_p(1,2*pi*R(i)*sqrt(eps_host(j))/lambda(j))-sqrt(eps_host(j))*riccati_bessel_psi(1,2*pi*R(i)*sqrt(eps_host(j))/lambda(j))*riccati_bessel_psi_p(1,2*pi*R(i)*sqrt(eps_filler(i,j))/lambda(j)))/(sqrt(eps_filler(i,j))*riccati_bessel_psi(1,2*pi*R(i)*sqrt(eps_filler(i,j))/lambda(j))*riccati_bessel_xi_p(1,2*pi*R(i)*sqrt(eps_host(j))/lambda(j))-sqrt(eps_host(j))*riccati_bessel_xi(1,2*pi*R(i)*sqrt(eps_host(j))/lambda(j))*riccati_bessel_psi_p(1,2*pi*R(i)*sqrt(eps_filler(i,j))/lambda(j)));
        RHS_i(i,j) = f(i)/R(i)^3*3i*lambda(j)^3/(16*pi^3*eps_host(j)^(3/2))*alpha_1i(i,j);
    end
end
RHS = sum(RHS_i,1);
epsilon_eff = transpose(eps_host).*(1+2*RHS)./(1-RHS);
end