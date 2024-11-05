%% function h_nu_z = struveh(nu,z)
% Simple implementation of the Struve function according to the integral
% representation (Only valid for Real nu > -1/2)
%
% H_{nu}(z) = 2*(z/2)^nu / ( sqrt(pi)*gamma(nu+1/2) ) *...
%               integral{from 0, to 1}( dt, (1-t.^2).^(nu-1/2) .* sin(z*t) )
%
% Ref: https://dlmf.nist.gov/11
%
% Inputs
%   * nu: (Scalar) Complex order, ( NOTE: Definition used is valid only for Re(nu) > -1/2 ) 
%   * z: (Vector or vector transpose) Complex variable
% 
% Outputs: 
% h_nu_z: (Real) vector
%
% Ex: 
% x = -10:0.01:10;
% n = 0:5; 
% y = zeros(length(n),length(x));
% for cont = 1:length(n)
%     y(cont,:) = struveh(n(cont),x);
% end
% plot(x,y,'linewidth',2); grid on;
% title('H_{nu}(z)');
% legend(num2str(n'));
% 
% Ex:
% x = 0:0.01:300;
% plot(x,struveh(0,x,'series'), x,struveh(0,x,'integral'), x,struveh(0,x,'integral2'))
% plot(x,struveh(1,x,'series'), x,struveh(1,x,'integral'), x,struveh(1,x,'integral2'))
%
% Ex: 
% x = -4:0.1:4;
% y = -4:0.1:4;
% [xx, yy] = meshgrid(x,y);
% zz = complex(xx,yy);
% H_0_z = struveh(0,zz);
% H_0_z = reshape(H_0_z,size(zz));
% figure, subplot 231
% surf(xx,yy,real(H_0_z))
% subplot 232
% surf(xx,yy,imag(H_0_z))
% subplot 233
% surf(xx,yy,abs(H_0_z))
% subplot 234
% imagesc(x,y,real(H_0_z))
% subplot 235
% imagesc(x,y,imag(H_0_z))
% subplot 236
% imagesc(x,y,abs(H_0_z))
%
% Germán A. Ramírez A.
% EPFL, MAG, August 2024

function H_nu_z = struveh(nu,z,varargin)
    definition = 'series';
    if ~isempty(varargin)
        definition = varargin{1};
    end

	z = z(:).';                     % Ensure proper vector transpose for multiplication with vector t
    
    if strcmpi(definition,'integral')
        t = linspace(0,0.999,1000)';    % Compute-intensive calculation. Do not use for large 'z' matrices! 

        % This integral expression is not advised due to oscillatory behavior
        H_nu_z = 2*(z/2).^nu ./ ( sqrt(pi).*gamma(nu+1/2) ) .*...
            trapz(t, (1-t.^2).^(nu-1/2) .* sin(t*z) );    

    elseif strcmpi(definition,'integral2')
        t1 = 0:0.001:0.1;
        t2 = 0.11:0.01:1;
        t3 = 1.1:0.1:10;
        t4 = 11:1:100;
        t = [t1,t2,t3,t4]';
        if nu == 0
            K_nu_z = 2/pi*trapz(t, exp(-sinh(t)*z));
        else
            K_nu_z = 2*(z/2).^nu ./ ( sqrt(pi).*gamma(nu+1/2) ) .*...
                trapz(t, exp(-t*z).*(1+t.^2).^(nu-1/2) );
        end
        Y_nu_z = bessely(nu,z);
        H_nu_z = K_nu_z + Y_nu_z;

    elseif strcmpi(definition,'series') 
        H_nu_z = zeros(1,length(z));
        zL = 20;       
        if ~isempty(z(z<zL))
            n = (0:50)';                     % Auto sum maximum (residual-based) should be implemented...  
            s = (-1).^n.*(1/2*z(z<zL)).^(2*n)./(gamma(n+3/2).*gamma(n+nu+3/2));
            H_nu_z(z<zL) = (1/2*z(z<zL)).^(nu+1).*sum(s,1); % Sum along dimension of n
        end
        if ~isempty(z(z>=zL))
            n = (0:10)';                     % Auto sum maximum (residual-based) should be implemented...  
            s = gamma(n+1/2).*(z(z>=zL)/2).^(nu-2*n-1)./gamma(nu+1/2-n);
            H_nu_z(z>=zL) = (1/pi).*sum(s,1) + bessely(nu,z(z>=zL)); 
        end
    end
end