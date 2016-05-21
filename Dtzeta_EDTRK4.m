function [nlin_z] = Dtzeta_EDTRK4(zhat)

global del2 KX KY;

% if abs(zhat(KX==0&KY==0))>1e-10
%     error('(kx,ky)=0');
% end

psihat = zhat./del2;

dpsidx  = real(ifft2(1i*KX.*psihat));
dpsidy  = real(ifft2(1i*KY.*psihat));
dzetadx = real(ifft2(1i*KX.*zhat));
dzetady = real(ifft2(1i*KY.*zhat));

J  = dpsidx.*dzetady - dpsidy.*dzetadx;
nlin_z = fft2( -J );

end