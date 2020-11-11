clear all;
tempo=clock;
% For simulation time
% calculation
%-----------------------------General Parameters---------------------------
Rb=40e9;
% Bit-Rate
Tb=1/Rb;
% Bit interval
samples=32;
% Samples per bit
Nb=64;
% Number of bits
t=Tb/samples*(0:samples-1);
time=Tb/samples*(0:Nb*samples-1);
f=1/(Nb*Tb)*([1:(samples*Nb)]-samples*Nb/2-1); % Simulation Bandwidth
omega=2*pi*f;
max_sections=20;
% Number of sections
BoptRX=1.35*Rb;
% Optical filter RX [Hz]
Be=0.8*Rb;
% Electrical filter [Hz]
%--------------------------------------------------------------------------
%-----------------------------Fiber Parameters-----------------------------
k=0.01;
% Stepsize constant
CR=1;
% Compesation ratio
% Standard Single Mode Faser (SSMF)
beta2 = -21.17e-27;
gama = 1.52e-3;
attenuation = (0.2*1e-3)/(10*log10(exp(1)));
L = 80e3;
%
%
%
%
Parameter for GVD [s^2/m]
Nonlinear parameter[1/(W*m)]
Attenuation [dB/m]
Fiber Length [m]
% Dispersion Compensating Fiber (DCF)
beta2_dcf = 131.88e-27;
% GVD Parameter [s^2/m]
gama_dcf = 5.27e-3;
% Nonlinear parameter[1/(W*m)]
attenuation_dcf = (0.50*1e-3)/(10*log10(exp(1))); % Attenuation [dB/m]
L_dcf = CR * L * abs(beta2) / abs(beta2_dcf);
% fiber length [m]
%--------------------------------------------------------------------------
%-----------------------------Transmitter----------------------------------
% Pulse 50% RZ
Pp=0.001;
% Pulse power [W]
p=sqrt(2*Pp)*cos(cos(pi*t/(Tb)).^2*pi/2);
% Precoding
a0=[0 0 0 0 0 0 1 0 0 0 0 1 1 0 0 0 1 0 1 0 0 1 1 1 1 0 1 0 0 0 1 1 ...
1 0 0 1 0 0 1 0 1 1 0 1 1 1 0 1 1 0 0 1 1 0 1 0 1 0 1 1 1 1 1 1];
% Sequence of pulses
u0=zeros(1,Nb*samples);
for j=1:Nb
u0( (j-1)*samples+1:j*samples )=p*a0(j);
end

P0=u0*u0’/length(u0);
%--------------------------------------------------------------------------
%-----------------------------General System Parameters--------------------
%-----------------------------Post-Compensation--------------------------
Pssmf_pos=10.^(-0.05)*1e-3;
% Optimum power for 40Gbit/s
stepsize_pos = L/ceil(L*gama*2*Pssmf_pos/k);% Stepesize in [m]
index=find(stepsize_pos>L);
% The stepsize should not
stepsize_pos(index)=L;
%exceed the fiber length
Pdcf_pos=10.^(-0.7)*1e-3;
% Optimum power for 40Gbit/s
stepsize_dcf_pos = L_dcf/ ...
ceil(L_dcf*gama_dcf*2*Pdcf_pos/k);
% Stepesize in [m]
index=find(stepsize_dcf_pos>L_dcf);
% The stepsize should not
stepsize_dcf_pos(index)=L_dcf;
%exceed the fiber length
u_pos=u0*sqrt(Pssmf_pos/P0);
% Input power
%--------------------------------------------------------------------------
%-----------------------------Pre-Compensation--------------------------
Pssmf_pre=10.^(-0.15)*1e-3;
% Optimum power for 40Gbit/s
stepsize_pre = L/ceil(L*gama*2*Pssmf_pre/k);% Stepesize in [m]
index=find(stepsize_pre>L);
% The stepsize should not
stepsize_pre(index)=L;
%exceed the fiber length
Pdcf_pre=10.^(-0.55)*1e-3;
% Optimum power for 40Gbit/s
stepsize_dcf_pre = L_dcf/ ...
ceil(L_dcf*gama_dcf*2*Pdcf_pre/k);
% Stepesize in [m]
index=find(stepsize_dcf_pre>L_dcf);
% The stepsize should not
stepsize_dcf_pre(index)=L_dcf;
%exceed the fiber length
u_pre=u0*sqrt(Pssmf_pre/(exp(attenuation*L)*P0)); % Input power
%--------------------------------------------------------------------------
%-----------------------------Hybrid-Compensation--------------------------
Pssmf_hyb=10.^(-0.1)*1e-3;
% Optimum power for 40Gbit/s
stepsize_hyb = L/ceil(L*gama*2*Pssmf_hyb/k);% Stepesize in [m]
index=find(stepsize_hyb>L);
% The stepsize should not
stepsize_hyb(index)=L;
%exceed the fiber length
Pdcf_hyb=10.^(-0.6)*1e-3;
% Optimum power for 40Gbit/s
stepsize_dcf_hyb = L_dcf/ ...
ceil(L_dcf*gama_dcf*2*Pdcf_hyb/k);
% Stepesize in [m]
index=find(stepsize_dcf_hyb>L_dcf);
% The stepsize should not
stepsize_dcf_hyb(index)=L_dcf;
%exceed the fiber length
u_hyb=u0*sqrt(Pdcf_hyb/(exp(attenuation_dcf*L_dcf/2)*P0)); % Input power
%--------------------------------------------------------------------------
U_pos=[u_pos];
U_pre=[u_pre];
U_hyb=[u_hyb];

% Post-Compensation--------------------------------------------------------
j=1;
while j<=max_sections
disp(sprintf(’Post-Comp. Span: %g Pdcf: %g dBm Pssmf: %g dBm’,...
j,10*log10(Pdcf_pos/1e-3),(10*log10(Pssmf_pos/1e-3))));
% Fiber SSMF
u_pos=splitstep(u_pos,L,stepsize_pos,attenuation,beta2,omega,gama);
% Optical Amplifier
Gain=Pdcf_pos*exp(attenuation*L)/Pssmf_pos;
u_pos=u_pos*sqrt(Gain);
% Fiber DCF
u_pos=splitstep(u_pos,L_dcf,stepsize_dcf_pos,attenuation_dcf, ...
beta2_dcf,omega,gama_dcf);
% Optical Amplifier
Gain=Pssmf_pos*exp(attenuation_dcf*L_dcf)/Pdcf_pos;
u_pos=u_pos*sqrt(Gain);
U_pos=[U_pos;u_pos];
j=j+1;
end
%--------------------------------------------------------------------------
% Pre-Compensation---------------------------------------------------------
j=1;
while j<=max_sections
disp(sprintf(’Pre-Comp. Span: %g Pdcf: %g dBm Pssmf: %g dBm’, ...
j,10*log10(Pdcf_pre/1e-3),(10*log10(Pssmf_pre/1e-3))));
% Optical Amplifier
Gain=Pdcf_pre*exp(attenuation*L)/Pssmf_pre;
u_pre=u_pre*sqrt(Gain);
% Fiber DCF
u_pre=splitstep(u_pre,L_dcf,stepsize_dcf_pre,attenuation_dcf, ...
beta2_dcf,omega,gama_dcf);
% Optical Amplifier
Gain=Pssmf_pre*exp(attenuation_dcf*L_dcf)/Pdcf_pre;
u_pre=u_pre*sqrt(Gain);
% Fiber SSMF
u_pre=splitstep(u_pre,L,stepsize_pre,attenuation,beta2,omega,gama);
U_pre=[U_pre;u_pre];
j=j+1;
end

%--------------------------------------------------------------------------
% Hybrid-Compensation------------------------------------------------------
j=1;
while j<=max_sections
disp(sprintf(’Hybrid-C. Span: %g Pdcf: %g dBm Pssmf: %g dBm’, ...
j,10*log10(Pdcf_hyb/1e-3),(10*log10(Pssmf_hyb/1e-3))));
% Fiber DCF/2
u_hyb=splitstep(u_hyb,L_dcf/2,stepsize_dcf_hyb/2, ...
attenuation_dcf,beta2_dcf,omega,gama_dcf);
% Optical Amplifier
Gain=Pssmf_hyb*exp(attenuation_dcf*L_dcf)/Pdcf_hyb;
u_hyb=u_hyb*sqrt(Gain);
% Fiber SSMF
u_hyb=splitstep(u_hyb,L,stepsize_hyb,attenuation,beta2,omega,gama);
% Optical Amplifier
Gain=Pdcf_hyb*exp(attenuation*L)/Pssmf_hyb;
u_hyb=u_hyb*sqrt(Gain);
% Fiber DCF/2
u_hyb=splitstep(u_hyb,L_dcf/2,stepsize_dcf_hyb/2, ...
attenuation_dcf,beta2_dcf,omega,gama_dcf);
U_hyb=[U_hyb;u_hyb];
j=j+1;
end
%--------------------------------------------------------------------------
MrkSz=16;
LnWdth=2;
plot(0:max_sections,eyeopeningpen(U_hyb,a0,Nb,samples,f,BoptRX,Be),...
’-o’,’LineWidth’,LnWdth,’MarkerSize’,MrkSz)
hold on
plot(0:max_sections,eyeopeningpen(U_pos,a0,Nb,samples,f,BoptRX,Be),...
’:sr’,’LineWidth’,LnWdth,’MarkerSize’,MrkSz)
plot(0:max_sections,eyeopeningpen(U_pre,a0,Nb,samples,f,BoptRX,Be),...
’--*k’,’LineWidth’,LnWdth,’MarkerSize’,MrkSz)
legend(’Hybrid Compensation’,’Post-Compensation’,’Pre-Compensation’,2);
set(gca,’linewidth’,3,’fontsize’,26);
xlabel(’Sections’,’FontSize’,30);
set(gca,’XTick’,0:2:max_sections)
ylabel(’Eye Opening Penalty [dB]’,’FontSize’,30);
title(sprintf(’RZ-OOK %gGbit/s’,Rb/1e9),’FontSize’,30);
grid on;
etime(clock,tempo)