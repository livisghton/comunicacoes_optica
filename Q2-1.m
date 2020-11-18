#Resoler a equação não linear de Schrodinger usando o Split-Step Fourier Method

#Ain e Aout são as amplitudes normalizadas da entrada e da saida
#L é o comprimento da fibra
#attenuation é "Alpha"
#beta2 é "Beta 2"
#omega é a frequencia angular
#gama é "Gama"
#P0 é a potencia do pulso
#LD é o comprimento de dispersão
#beta2 é a dispersão anomala
#Tq é uma medida da largura de pulso

#P0 * LD * gama = 1

beta3 = 0;

LD = Tq^2/mod(beta2)

function res = A(z,t)
    res = sqrt(P0 * sech(t/Tq)*exp((z/2*LD)*(-j))) 
endfunction


function Aout = splitstep(Ain, L, stepsize, attenuation, beta2, omega, gama)