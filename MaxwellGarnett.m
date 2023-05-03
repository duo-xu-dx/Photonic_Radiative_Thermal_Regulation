function epsilon_eff = MaxwellGarnett(lambda,num_filler,f,eps_host,eps_filler)
%MAXWELLGARNETT Computes the effective dielectric constant using Maxwell-Garnett theory.
%   epsilon_eff = MAXWELLGARNETT(lambda, num_filler, f, eps_host, eps_filler) calculates the effective 
%   dielectric constant (epsilon_eff) for a composite material using the Maxwell-Garnett effective medium 
%   approximation. The function takes the following inputs:
%   
%   lambda:     Wavelengths as a vector (not used in the calculations, but required for consistency with 
%               the MaxwellGarnettMie function).
%   num_filler: Number of filler materials in the composite.
%   f:          Volume fractions of the filler materials as a vector of length num_filler.
%   eps_host:   Dielectric constant of the host material as a scalar.
%   eps_filler: Dielectric constants of the filler materials as a vector of length num_filler.
%
%   epsilon_eff: The effective dielectric constant of the composite material as a vector with the same length as lambda.
%
%   The function calculates the effective dielectric constant using the Maxwell-Garnett effective medium
%   approximation, which models the composite material as a host material with embedded spherical
%   inclusions (filler materials).

RHS_i = zeros(num_filler,numel(lambda));
for i = 1:num_filler
    for j = 1:numel(lambda)
        RHS_i(i,j) = f(i)*(eps_filler(i,j)-eps_host(j))/(eps_filler(i,j)+2*eps_host(j));
    end
end
RHS = sum(RHS_i,1);
epsilon_eff = transpose(eps_host).*(1+2*RHS)./(1-RHS);
end