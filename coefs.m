function [a0,an,bn] = coefs(y)
% returns fourier coefficients, input column vector
nx=length(y);
xs=fftshift(y);
fspace=(2/nx)*fft(xs);
a0=fspace(1);
an=real(fspace(2:nx/2+1));
bn=-imag(fspace(2:nx/2+1));
end
