function y = eyeopeningpen(U_hyb,a0,Nb,samples,f,BoptRX,Be)
    GF = gaussian(f,BoptRX);
    GF_F = fftshift(GF);
    BF = besselMC(f,Be);
    BF_F = fftshift(BF);
    [line,row] = size(U_hyb);
    Eop=[];
    for i=1:line
        U = U_hyb(i,:); % Linha avaliada
        UF = fft(U);    % Transformação para dominio da Frequencia
        U_GF = UF.*GF_F;  % Filtro Gaussiano
        U = ifft(U_GF);
        U_detected = abs(U).^2; %Detção do Sinal
        U_detected_F = fft(U_detected);
        U_BF = U_detected_F.*BF_F; %Filtro de Bessel
        U_T = ifft(U_BF); %Transformação dominio do tempo
        U_Final = real(U_T);
        
        Min1 = realmax;
        Max0 = realmin;
        centralSample = ceil(samples/2);
        k=centralSample;
        for j=1:Nb
            
            if a0(j)==0
                if U_Final(k)>Max0
                    Max0 = U_Final(k);
                end
            else
                if U_Final(k)< Min1
                    Min1 = U_Final(k);
                end
            end
            k = k+samples;
            j=j+1;
        end
        Eop(i) = Min1-Max0;
    end
    for i=1:line
        y(i) = 10*log10(Eop(1)/Eop(i));
    end
end