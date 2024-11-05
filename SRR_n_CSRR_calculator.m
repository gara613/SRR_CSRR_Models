%% function [f0_SRR, L, C] = SRR_n_CSRR_calculator(rext, c, d, eps_r, h_subs)
%
% Juan Domingo Baena Doello et al, "Equivalent-Circuit Models for Split-Ring Resonators 
% and Complementary Split-Ring Resonators Coupled to Planar Transmission Lines"
% IEEE TRANSACTIONS ON MICROWAVE THEORY AND TECHNIQUES, VOL. 53, NO. 4, APRIL 2005 
%
% R. Marqués, F. Mesa, J. Martel, and F. Medina, “Comparative analysis of edge- 
% and broadside-coupled split ring resonators for metamaterial design—Theory and experiment,” 
% IEEE Trans. Antennas Propag., vol.51, no. 10, pp. 2572–2581, Oct. 2003. 
%
% SRR: EC (Edge coupled)
% CSRR: 
% 
% * Inputs:  
%   - rext:     (scalar) External radius
%   - c:        (scalar) Trace or Slot width (SRR / C-SRR)
%   - d:        (scalar) Distance between traces/slots
%   - eps_r:    (vector) Relative permittivity of substrate
%   - h_subs:   (vector) Substrate height
%
% * Outputs:
%   - f0_SRR:   (Struct) with fields 'f0_EC_SRR' and 'f0_C_SRR' each is a matrix of size [length(eps_r),length(h_subs)]
%   - L:        (Struct) with fields 'Ls' and 'Lc', each representing the inductance of the SRR/CSRR  
%   - C:        (Struct) with fields 'Cs' and 'Cc' each is a matrix of size [length(eps_r),length(h_subs)] representing the capacitance of the SRR/CSRR  
%
% Ex: delta(f_r) < 1% w.r.t Juan's code (written in Fortran).
% rext = 2.3e-3;
% c = 0.2e-3;
% d = 0.2e-3;
% eps_r = 2.3; 
% h_subs = 0.5e-3; 
% [f0_SRR, L, C] = SRR_n_CSRR_calculator(rext, c, d, eps_r, h_subs)
% 
% Ex: Recreate Fig. 6 from Juan Domingo's paper 
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
% Germán A. Ramírez 
% MAG - EPFL, July 2024

function [f0_SRR, L, C] = SRR_n_CSRR_calculator(rext, c, d, eps_r, h_subs, varargin)
    mu0 = 4*pi*1e-7;
    eps0 = 8.8541878188*1e-12;
    c0 = 1/sqrt(mu0*eps0);
    
    % Varargin: (Experimental) Superstrate permittivity. Accuracy not tested. 
    % Use instead "SRR_n_CSRR_ML_calc" for multilayer structures.   
    if ~isempty(varargin)
        epsr_ss = varargin{1};
    else
        epsr_ss = 1; 
    end

    r0 = rext - c - d/2;

    k1 = 0.01:0.01:1;
    k2 = 1.1:0.1:10;
    k3 = 11:1:100;
    k4 = 110:10:10000;
    k5 = 10100:100:1e5;
    k = [k1,k2,k3,k4,k5];

    % Inductance of SRR
    a = r0 - c/2; 
    b = r0 + c/2; 
    integ = trapz(k, 1./k.^2.*( b.*B(k*b) - a.*B(k*a) ).^2);
    L.Ls = mu0*pi^3/(4*c^2)*integ;

    % Inductance of CSRR
    aa = d;
    bb = d+2*c;
    kk = aa/bb; 
    L_pul = mu0/(4*K_Kp(kk));
    L0 = 2*pi*r0*L_pul;
    L.Lc = L0/4;
    
    f0_SRR.f0_EC_SRR = zeros(length(eps_r),length(h_subs));
    f0_SRR.f0_C_SRR = zeros(length(eps_r),length(h_subs));
    C.Cs = zeros(length(eps_r),length(h_subs));
    C.Cc = zeros(length(eps_r),length(h_subs));

    for cont_eps_r = 1:length(eps_r)
        for cont_h_subs = 1:length(h_subs)
            %% SRR: R. Marqués
            [Z,eps_eff] = Z0_CPS(d,c,h_subs(cont_h_subs),eps_r(cont_eps_r));

            C_pul = sqrt(eps_eff)/(c0*Z);
            C0 = 2*pi*r0*C_pul;
            C.Cs(cont_eps_r,cont_h_subs) = C0/4;
            f0_SRR.f0_EC_SRR(cont_eps_r,cont_h_subs) = sqrt(1/(L.Ls*C.Cs(cont_eps_r,cont_h_subs)))/(2*pi);

            %% C-SRR: J.D.B.D
            integSubs = trapz(k, (1./k.^2.*( b.*B(k*b) - a.*B(k*a) ).^2).*...
                (1/2*(1 + (1+eps_r(cont_eps_r)/epsr_ss*tanh(k*h_subs(cont_h_subs))) ./ (1+epsr_ss/eps_r(cont_eps_r)*tanh(k*h_subs(cont_h_subs))) )) );

            C.Cc(cont_eps_r,cont_h_subs) = eps0*pi^3/(c^2)*integSubs;           
            f0_SRR.f0_C_SRR(cont_eps_r,cont_h_subs) = sqrt(1/(L.Lc*C.Cc(cont_eps_r,cont_h_subs)))/(2*pi);
        end
    end
end


%%
function y = B(x)
    y = struveh(0,x,'series').*besselj(1,x) - struveh(1,x,'series').*besselj(0,x);
end

%% Stripline impedance from Microwave Solid State Circuit Design, 2nd Edition, Chapter 2
function y = K_Kp(k)
    kp = sqrt(1-k.^2);
    y = zeros(length(k),1);
    
    y(k>=0 & k<=1/sqrt(2)) = pi./log( 2*(1+sqrt(kp))./(1-sqrt(kp)) ) ;
    y(k>=1/sqrt(2) & k<=1) = 1/pi*log( 2*(1+sqrt(k))./(1-sqrt(k)) ) ;
end

% W-S-W: CPW and CPS as complementary structures
% For the CPW 
% W: separation between traces
% S: Central Trace width
function [Z,eps_eff] = Z0_CPW(S,W,h,eps_r)
    a = S/2;
    b = S/2+W;
    k = a/b; 
    k1 = sinh(pi*a./(2*h))./sinh(pi*b./(2*h));
    
    eps_eff = 1+(eps_r-1)/2*K_Kp(k1)./K_Kp(k);    
    Z = 30*pi/sqrt(eps_eff)./K_Kp(k);
end

% W-S-W: CPW and CPS as complementary structures
% For the CPS
% W: Traces width
% S: Separation between traces
function [Z,eps_eff] = Z0_CPS(S,W,h,eps_r)
    a = S/2;
    b = S/2+W;
    k = a/b; 
    k1 = sinh(pi*a./(2*h))./sinh(pi*b./(2*h));
    
    eps_eff = 1+(eps_r-1)/2*K_Kp(k1)./K_Kp(k);    
    Z = 120*pi/sqrt(eps_eff).*K_Kp(k);
end