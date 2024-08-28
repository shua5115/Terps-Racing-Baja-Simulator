function [Tp1,Tp0,Cpx,Fpx,Fps,Fflyx,Ffly,wpc,MPMax] = ...
    PrimaryV3(n,Tp,MBelt,Ts1,Wp,X1,Rp,Alpha,BeltSH,Uep,Ufly,Theta,Kp,Xp0,Delta,Mfly,Rfly,Mlink,Rlink,A,B,C)
%Calculate Primary RPM
%Taught side tenion
Tp1=Ts1;
%Slack side tenion
Tp0=Tp1./exp(Uep.*Alpha);
%Accounting for Motor Torque
Tp1=Tp1+((Tp.*6)./Rp);
Tp0=Tp0-((Tp.*6)./Rp);
%Tension Component in X Direction
Tp0x=zeros(1,n);
for i=1:1:n
    if Alpha(i)<=pi
        Tp0x(i)=Tp0(i).*cos(.5*(pi-Alpha(i)));
    elseif Alpha(i)>pi
        Tp0x(i)=Tp0(i).*cos(.5*(Alpha(i)-pi));
    end
end
Tp1x=zeros(1,n);
for i=1:1:n
    if Alpha(i)<=pi
        Tp1x(i)=Tp1(i).*cos(.5*(pi-Alpha(i)));
    elseif Alpha(i)>pi
        Tp1x(i)=Tp1(i).*cos(.5*(Alpha(i)-pi));
    end
end
%Centrifugal Force
Wp=Wp*(pi/30);
Cp=(1/12).*MBelt.*(Rp.^2).*(Wp.^2);
%Distributed Centrifugal Force
Cpx=2.*Cp.*sin(Alpha./2); %Total Centrifugal Force
%Required side force
Fpx=(Tp0x+Tp1x-Cpx); %Total clamping force
Fp=(Fpx./(2.*sin(Alpha./2)))./(2.*tand(Theta/2)); %Distributed clamping force
%Pressure Spring
Fps=Kp(A).*(Xp0+X1); %Pressure spring force vector (lbf)
%Required Flyweight force
Fflyx=Fp+Fps;
Ffly=(Fflyx./tand(Delta(B,:)))-(Fflyx./tan(acos((Rfly-1.625)./1.25)));
%Required engine speed
wpc=(12.*Ffly./((Mfly(C).*Rfly)+(Mlink.*Rlink))).^.5;
wpc=(wpc.*30)./(pi); %Rev/Min
%Maximum Torque Tranferable (Without Slip)
%MPMax=(Tp1-Tp0).*(Rp./12);
%MPMax=(Tp1xTp0x)./cosd(Theta/2).*Uep.*(Rp./12);
MPMax=(Fflyx./cosd(Theta/2)).*Uep.*(Rp./12);
%Ft*lbs
end
