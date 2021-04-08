%-------------------------------------------------------------------------------
% Satillite Simulation of Pitch Attitude Maneuver - P.J.Barnhart/CWRU - 2015
%-------------------------------------------------------------------------------
function satsim = satsim(t,gains) % input time vector and controller gains
%
wn=0.23*2*pi; wn2=wn*wn; zeta=0.569; eta=3.445; Ist=275; Irw=0.995; % constants
%
kp=gains(1)*(wn2/Ist); % controller proportional gain
ki=gains(2)*(wn2/Ist); % controller integral gain
kd=gains(3)*(wn2/Ist); % controller derivative gain
%
num=[0 0 0 kd kp ki]; % Closed-Loop Transfer Function
den=[1 2*zeta*wn (wn2+kd*eta) (kp*eta+kd) (ki*eta+kp) ki];
%
figure(1); % Pitch Attitude
a=12; % 12 degree step command for attitude maneuver
c=step(a*num,den,t);
[mt,nt]=size(t); 
iset=0; 
imax=0; 
cmax=0.0; % find settling time & overshoot
for ni = 1:nt
    if (abs(c(ni)-a) > 0.2) 
        iset=ni; 
    end
    if (c(ni) > cmax) 
        imax=ni; 
        cmax=c(ni); 
    end
end
if (iset < nt)
    iset=iset+1; 
end
settling=t(iset); 
overshoot=c(imax)/a-1; 
ymin=0; 
ymax=18;
plot(t,c,'b', [t(1) t(nt)],[a a],'m:',[t(1) t(nt)],[1.2*a 1.2*a],'r-', [t(iset) t(iset)],[ymin ymax],'g-');
axis([0.0 t(mt,nt) ymin ymax]);
xlabel('Time (s)'); ylabel('Pitch Attitude (deg)');
%
figure(2); % Closed-Loop Bode Plot
[mag,phs,w]=bode(num,den);
[mb,nb]=size(mag); ibw=0; % find bandwidth
for ni = 1:mb
   if (abs(mag(ni)) > 0.707)
       ibw=ni; 
   end
end
bw=w(ibw+1);
subplot(2,1,1); semilogx(w,20*log10(mag));
xlabel('Frequency (rad/sec)'); ylabel('Gain (dB)'); grid on;
subplot(2,1,2); semilogx(w,phs);
xlabel('Frequency (rad/sec)'); ylabel('Phase (deg)'); grid on;
%
figure(3); % Root-Locus Plot
rlocus(num,den);
%
num2=[0 0 kd kp ki 0]; % Angular Rate
cdot=step(a*num2,den,t);
%
figure(4); % Satellite Momentum
plot(t,Ist*cdot*pi/180,'b',[t(1) t(nt)],[0 0],'m:', ...
     [t(1) t(nt)],[52.10 52.10],'r-', ...
     [t(1) t(nt)],[-52.10 -52.10],'r-');
xlabel('Time (s)'); ylabel('Satellite Momentum (lb.ft.s)');
%
figure(5); % Reaction Wheel Speed
plot(t,(Ist/Irw)*cdot*(pi/180)*(30/pi),'b',[t(1) t(nt)],[0 0],'m:', ...
     [t(1) t(nt)],[500 500],'r-', ...
     [t(1) t(nt)],[-500 -500],'r-');
xlabel('Time (s)'); ylabel('Reaction Wheel Speed (rpm)');
%
num3=[0 kd kp ki 0 0]; % Angular Acceleration
cdbldot=step(a*num3,den,t);
%
figure(6); % Reaction Wheel Motor Torque
Tqmax=6.565;
plot(t,Ist*cdbldot*(pi/180),'b',[t(1) t(nt)],[0 0],'m:', ...
     [t(1) t(nt)],[Tqmax Tqmax],'r-', ...
     [t(1) t(nt)],[-Tqmax -Tqmax],'r-');
xlabel('Time (s)'); ylabel('Actuator Torque (lb.ft)');
%
figure(7); % Total Kinetic Energy
plot(t,0.5*Ist*(cdot.*cdot)*(pi/180)*(pi/180)*(1+(Ist/Irw)),'b')
xlabel('Time (s)'); ylabel('Kinetic Energy, Satellite & RW (lb.ft)');
% find maximum values for RWA control authority
icd=0; icdd=0; icddd=0; icdm=0; icddm=0; icdddm=0; cdm=0.0; cddm=0.0; cdddm=0.0;
for ni = 1:nt
   if (abs(cdot(ni)) > cdm) icdm=ni; cdm=cdot(ni); end
   if (abs(cdbldot(ni)) > cddm) icddm=ni; cddm=cdbldot(ni); end
   if (abs(cdot(ni)*cdbldot(ni)) > cdddm)
      icdddm=ni; cdddm=cdot(ni)*cdbldot(ni);
   end
end
maxWw=(Ist/Irw)*cdot(icdm)*(pi/180)*(30/pi);
maxTq=Ist*cdbldot(icddm)*(pi/180);
maxHp=(Ist)*cdbldot(icdddm)*(Ist/Irw)*cdot(icdddm)*(pi/180)*(pi/180)/550;
%
numol=[0 0 kd*eta (kp*eta+kd) (ki*eta+kp) ki]; % Open-Loop TF
denol=[1 2*zeta*wn wn2 0 0 0];
%
figure(8); % Open-Loop Bode Plot
[mag,phs,w]=bode(numol,denol);
subplot(2,1,1); semilogx(w,20*log10(mag));
xlabel('Frequency (rad/sec)'); ylabel('Gain (dB)'); grid on;
subplot(2,1,2); semilogx(w,phs);
xlabel('Frequency (rad/sec)'); ylabel('Phase (deg)'); grid on;
% Integrations to find ITAE and Total Energy
[m,n]=size(c');
[aa]=a*ones(m,n);
[y]=[t].*abs([c']-[aa]);
ITAE=trapz(t,y);
[z]=Ist*([cdot'].*[cdbldot'])*(1+(Ist/Irw))*(pi/180)*(pi/180);
for ni = 1:nt
   if (ni > iset) z(ni)=0; end
end
Etot=trapz(t,z);
%
disp('    Overshoot Settling  Speed     Torque    Hp        BW        Energy');
disp([overshoot settling maxWw maxTq maxHp bw Etot]);
satsim=ITAE;
return
end
