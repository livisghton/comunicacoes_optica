clear all;
tempo=clock;
% For simulation time
% calculation
%-----------------------------General Parameters---------------------------
Rb=10e9;
% Bit rate [bit/s]
Tb=1/Rb;
% Bit interval [s]
q0=3.5;
% Separation between solitons
Tq=1/(2*q0*Rb);
% Soliton width [s]
comp=4;
% Number of bit slots
samples=64;
% Samples per slot
N=comp*samples;
% Total number of samples
time=Tb/samples*(-N/2:(N/2-1));
% Time vector for all slots
f=1/(comp*Tb)*( [1:N] - N/2 - 1);
% Simulation Bandwidth [Hz]
omega=2*pi*f;
% Angular frequency
%--------------------------------------------------------------------------
%-----------------------------Fiber Parameters-----------------------------
k=0.01;
% Stepsize constant
beta2 = -1e-27;
% Parameter for GVD [s^2/m]
gama = 1.52e-3;
% Nonlinear parameter[1/(W*m)]
attenuation = (0.0*1e-3)/(10*log10(exp(1))); % Attenuation [dB/m]
Ld=Tq^2/abs(beta2);
% Dispersion length [m]
L = 50*Ld;
% Fiber Length [m]
L_step=Ld;
% Step for the 3D-plot [m]
P0=abs(beta2)/(gama*Tq^2);
% Pulse peak power [W]
stepsize = L/ceil(L*gama*P0/k);
% Stepesize in [m]
index=find(stepsize>L);
% The stepsize should not
stepsize(index)=L;
%exceed the fiber length
%--------------------------------------------------------------------------
%-----------------------------Initial Condition----------------------------
% 0 - One soliton pulse
% 1 - One Gaussian pulse
% 2 - Soliton interaction (Two pulses)
% 3 - Soliton interaction (Two pulses with different amplitudes)
initial_pulse=3;
switch initial_pulse
case 0
u0=sqrt(P0)*sech(time/Tq);
case 1
u0=sqrt(P0)*exp(-time.^2/(2*Tq^2));
case 2
u0=sqrt(P0)*( sech((time-Tb/2)/Tq) + sech((time+Tb/2)/Tq));
case 3
u0=sqrt(P0)*( sech((time-Tb/2)/Tq) + 1.1*sech((time+Tb/2)/Tq));
end
%--------------------------------------------------------------------------
MrkSz=16;

LnWdth=2;
figure(1);
plot(time/1e-12,(abs(u0).^2)/1e-3,'LineWidth',LnWdth,'MarkerSize',MrkSz);
hold on;
% Propagating through the fiber
u=splitstep1(u0,L,stepsize,attenuation,beta2,omega,gama);
figure(1)
plot(time/1e-12,(abs(u).^2)/1e-3,'rx','LineWidth',LnWdth,'MarkerSize',MrkSz);
set(gca,'linewidth',3,'fontsize',26);
xlabel('Time [ps]','FontSize',30);
ylabel('Power [mW]','FontSize',30);
% 3D-plot
u=u0;
U=u0;
J=L/L_step;
for j=1:J
u=splitstep1(u,L_step,stepsize,attenuation,beta2,omega,gama);
U=[U;u];
end
figure(2);
h=waterfall(time/1e-12,(0:size(U,1)-1)*L_step/1e3,abs(U).^2/1e-3);
set(gca,'linewidth',3,'fontsize',26);
set(h,'LineWidth',LnWdth);
xlabel('Time [ps]','FontSize',30);
ylabel('Distance [Km]','FontSize',30);
zlabel('Power [mW]','FontSize',30);
axis([time(1)/1e-12 time(end)/1e-12 0 J*L_step/1e3 ...
0 max(max((abs(U).^2)./1e-3)) 15 16]);
set(gca,'YDir','reverse');
etime(clock,tempo)