function Aout = splitstep1(Ain, L, stepsize, attenuation, beta2, omega, gama)
    %beta2 = parametro de dispers�o de segunda ordem
    %beta3 = parametro de dispers�o de terceira ordem
    %alpha = attenuation = coeficiente de atenua��o na fibra
    %gama = coeficiente n�o linear
    %omega = frequencia angular
    %L = tamanho da fibra
    %Ain = Potencia de entrada
    
    j = sqrt(-1);
    omega = fftshift( omega );
    K = L / stepsize;%n�mero de divis�es
    for k = 1 : K
        n = -j * gama * (abs(Ain).^2); %operador n�o linear
        Ain = Ain.* exp(n * (stepsize / 2));
        uf = fft( Ain );
        d = -(attenuation / 2)-((j / 2) * beta2 * (omega.^2)); %operador linear
        uf = uf.* exp(d * stepsize);
        Ain = ifft(uf);
        n = -j * gama * (abs(Ain).^2);
        Ain = Ain.*exp(n * (stepsize / 2));
    end
    Aout = Ain;
end