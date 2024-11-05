%% function [f0_SRR, L, C] = SRR_n_CSRR_ML_calc(rext, c, d, eps_r, h_subs)
%
% Equivalent-Circuit Models for Split-Ring Resonators and Complementary Split-Ring Resonators 
% in a Multilayer Subsstrate stack-up
% Germán A. Ramírez, Anja K. Skrivervik, "Permittivity Sensors based on Conventional 
% and Complementary Split Ring Resonators", EuCAP 2025
%
% SRR: EC, BC
% CSRR: 
% 
% * Inputs:  
%   - rext:     (scalar) External radius
%   - c:        (scalar) Trace or Slot width (SRR / C-SRR)
%   - d:        (scalar) Distance between traces/slots
%   - eps_r:    (vector) Top-Bottom Relative permittivity of the substrates' stack-up, (without including the air layers assumed above and below)
%   - h_subs:   (vector) Top-Bottom Substrate height of the substrates' stack-up, (without including the air layers assumed above and below)
%
% * Outputs:
%   - f0_SRR:   (scalar) Struct with fields 'f0_EC_SRR' and 'f0_C_SRR' representing the resonant frequencies of the SRR and CSRR
%   - L:        (scalar) Struct with fields 'Ls' and 'Lc' representing the Inductance of the SRR and CSRR
%   - C:        (scalar) Struct with fields 'Cs' and 'Cc' representing the Capacitance of the SRR and CSRR
%
% Ex: Recreate Fig. 6. from Juan Domingo's paper 
% rext = 2.3e-3;
% c = 0.2e-3;
% d = 0.2e-3;
% eps_r = [2,4,6,8,10];
% h_subs = (0:0.01:2)*1e-3;
% f0_SRR = SRR_n_CSRR_calculator(rext, c, d, eps_r, h_subs);
% figure,
% plot(h_subs/1e-3,f0_SRR.f0_EC_SRR/1e9,'linewidth',2), grid on
% hold on,
% plot(h_subs/1e-3,f0_SRR.f0_C_SRR/1e9,'--','linewidth',2)
% mycolors = [ [0 0.4470 0.7410];
% [0.8500 0.3250 0.0980];
% [0.9290 0.6940 0.1250];
% [0.4940 0.1840 0.5560];
% [0.4660 0.6740 0.1880]; ];
% ax = gca; 
% ax.ColorOrder = mycolors;
% legend(num2str([eps_r,eps_r]'),'numColumns',2)
% xlabel('Dielectric thickness (mm)')
% ylabel('f_0 (GHz)')
% eps_r = 1:0.1:10;
% h_subs = (0.1:0.1:0.5)*1e-3;
% f0_SRR = SRR_n_CSRR_calculator(rext, c, d, eps_r, h_subs);
% figure,
% plot(eps_r,f0_SRR.f0_EC_SRR/1e9,'linewidth',2), grid on
% hold on,
% plot(eps_r,f0_SRR.f0_C_SRR/1e9,'--','linewidth',2)
% ax = gca; 
% ax.ColorOrder = mycolors;
% legend(num2str([h_subs*1e3,h_subs*1e3]'),'numColumns',2)
% xlabel('relative dielectric constant (\epsilon_r^'')')
% ylabel('f_0 (GHz)')
%
% Ex:
% See the script 'freqResonance_vs_MUT_eps_SRR.m' for generating the plots in the paper 
% Germán A. Ramírez, Anja K. Skrivervik, "Permittivity Sensors based on Conventional and Complementary Split Ring Resonators"
%
% Germán A. Ramírez 
% MAG - EPFL, September 2024

function [f0_SRR, L, C] = SRR_n_CSRR_ML_calc(rext, c, d, eps_r, h_subs, varargin)
    mu0 = 4*pi*1e-7;
    eps0 = 8.8541878188*1e-12;
    c0 = 1/sqrt(mu0*eps0);
    
    if ~isempty(varargin)
        switch length(varargin)
            case 1
                NU = varargin{1};           % Number of superstrates
            case 2
                NU = varargin{1};
                gnd_low = varargin{2};      % Define the presence of a lower ground plane
        end
    else
        NU = length(eps_r)-1;               % By defualt, the resonator is on top of first substrate
        gnd_low = 0;                        % No ground plane as part of the resonator
    end

    eps_r = [1; eps_r; 1];

    r0 = rext - c - d/2;

    k1 = 0.01:0.01:1;
    k2 = 1.1:0.1:10;
    k3 = 11:1:100;
    k4 = 110:10:10000;
    k5 = 10100:100:1e5;
    k = [k1,k2,k3,k4,k5];

    a = r0 - c/2; 
    b = r0 + c/2; 
    integ = trapz(k, 1./k.^2.*( b.*B(k*b) - a.*B(k*a) ).^2);
    L.Ls = mu0*pi^3/(4*c^2)*integ;

    %% SRR: 
    [Z,eps_eff] = Z0_CPS(d,c,h_subs,eps_r,NU);

    C_pul = sqrt(eps_eff)/(c0*Z);
    C0 = 2*pi*r0*C_pul;
    C.Cs = C0/4;
    f0_SRR.f0_EC_SRR = sqrt(1/(L.Ls*C.Cs))/(2*pi);

    %% C-SRR:  
    if ~gnd_low
        Yi = eps_r(2:end-1).* (1 + eps_r(2:end-1).*tanh(h_subs*k)) ./ (eps_r(2:end-1) + tanh(h_subs*k));
    else    % In normal cases Gnd plane is not part of the resonator and its presence should be accounted by the excitation mechanism
        Yi = eps_r(2:end-1).* (1 + eps_r(2:end-1).*tanh(h_subs*k)) ./ (eps_r(2:end-1) + tanh(h_subs*k));
        Yi(end,:) = eps_r(end-1).* coth(h_subs(end)*k);
    end
    
    % PENDING: General solution (Recursive use of Yi)
	% Not a very elegant implementation, but works for the three cases of interest...    
    N_diels = length(eps_r)-2;
    switch N_diels     
        case 1
            Y = 1/2*(1+Yi);
        case 2
            Y = 1/2*sum(Yi); 
        case 3
            Y = 1/2*...
            ( eps_r(3)* (Yi(1,:) + eps_r(3)*tanh(k*h_subs(2))) ./ (eps_r(3) + Yi(1,:).*tanh(k*h_subs(2))) +...
            Yi(3,:) );
    end  
        
    integSubs = trapz(k, (1./k.^2.*( b.*B(k*b) - a.*B(k*a) ).^2).*Y);

    C.Cc = eps0*pi^3/(c^2)*integSubs;
    aa = d;
    bb = d+2*c;
    kk = aa/bb; 
    L_pul = mu0/(4*K_Kp(kk));
    L0 = 2*pi*r0*L_pul;
    L.Lc = L0/4;
    f0_SRR.f0_C_SRR = sqrt(1/(L.Lc*C.Cc))/(2*pi);
end


%% Testing function (Same as used in Marques and Baena's works)
function y = B(x)
    y = struveh(0,x,'series').*besselj(1,x) - struveh(1,x,'series').*besselj(0,x);
end

%% Stripline impedance from "Microwave Solid State Circuit Design", 2nd Edition, Chapter 2
function y = K_Kp(k)
    kp = sqrt(1-k.^2);
    y = zeros(length(k),1);
    
    % Hilberg (Equivalent to the expression below)
%     y(k>=0 & k<=1/sqrt(2)) = 0.5*pi./log(2*sqrt(1+kp(k>=0 & k<=1/sqrt(2))) ./ sqrt(1-kp(k>=0 & k<=1/sqrt(2))) ) ;
%     y(k>=1/sqrt(2) & k<=1) = 2/pi*log( 2*sqrt( 1+k(k>=1/sqrt(2) & k<=1) )  ./ sqrt( 1-k(k>=1/sqrt(2) & k<=1) ) );
%     
    % Bahl Bhartia
	y(k>=0 & k<=1/sqrt(2)) = pi./log( 2*(1+sqrt(kp(k>=0 & k<=1/sqrt(2))))./(1-sqrt(kp(k>=0 & k<=1/sqrt(2)))) ) ;
    y(k>=1/sqrt(2) & k<=1) = 1/pi*log( 2*(1+sqrt(k(k>=1/sqrt(2) & k<=1)))./(1-sqrt(k(k>=1/sqrt(2) & k<=1))) ) ;
end

%% W-S-W: CPW and CPS as complementary structures
% For the CPS
% W: Traces width
% S: Separation between traces
function [Z,eps_eff] = Z0_CPS(S,W,h,eps_r,NU)
    a = S/2;
    b = S/2+W;
    k = sqrt(1-(a/b)^2); 
    k_i = sqrt( 1 - sinh(pi*a./(2*h)).^2 ./ sinh(pi*b./(2*h)).^2 );

    eps_eff = 1 + K_Kp(k)/2*sum( (eps_r(2:NU+1)-eps_r(1:NU))./K_Kp(k_i(1:NU)) )...
                + K_Kp(k)/2*sum( (eps_r(NU+2:end-1)-eps_r(NU+3:end))./K_Kp(k_i(NU+1:end)) );

    Z = 120*pi/sqrt(eps_eff)./K_Kp(k);
end

%% W-S-W: CPW and CPS as complementary structures
% For the CPW 
% W: separation between traces
% S: Central Trace width
function [Z,eps_eff] = Z0_CPW(S,W,h,eps_r)
    a = S/2;
    b = S/2+W;

	k = sqrt(1-(a/b)^2); 
    k_i = sqrt( 1 - sinh(pi*a./(2*h)).^2 ./ sinh(pi*b./(2*h)).^2 );

    eps_eff = 1 + K_Kp(k)/2*sum( (eps_r(2:NU+1)-eps_r(1:NU))./K_Kp(k_i(1:NU)) )...
                + K_Kp(k)/2*sum( (eps_r(NU+2:end-1)-eps_r(NU+3:end))./K_Kp(k_i(NU+1:end)) );
  
    Z = 30*pi/sqrt(eps_eff).*K_Kp(k);
end