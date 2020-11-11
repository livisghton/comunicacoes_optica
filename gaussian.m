function y=gaussian(f,Bopt)
% Returns the value of the 2nd order gaussian optical Filter
% for a given frequency and bandwidth

y=exp(-(f/(0.85*Bopt)).^2);
