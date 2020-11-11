function y=besselMC(f,Be)
% Returns a complex value of the 5th order bessel Electrical Filter
% for a given frequency and bandwidth. The constant delay is compensated by
% multiplying by "exp(i*2*pi*F.*tau)"

% Group Delay Compensation
% tau=1/(2*pi)*(893025 + 99225*(f*2*pi).^2 + 6300*(f*2*pi).^4 + ...
% 315*(f*2*pi).^6 + 15*(f*2*pi).^8)./(893025 + 99225*(f*2*pi).^2 + ...
% 6300*(f*2*pi).^4 + 315*(f*2*pi).^6 + 15*(f*2*pi).^8 + (f*2*pi).^10);
% f=0 
tau=1/(2*pi);

F=2.427410702*f/Be;
y=( 945./(i*F.^5 + 15*F.^4 - 105*i*F.^3 - 420*F.^2 ...
    + 945*i*F + 945) ).*exp(i*2*pi*F.*tau);
