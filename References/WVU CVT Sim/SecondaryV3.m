function [Fxts,Fyts,FHelix,Fs,Fsx,Csx,Ts0,Ts1,MSMax,Effs] = ...
    SecondaryV3(n,MBelt,Ws,Tv,Ts,X2,Y2,Rs,Phi,Beta,BeltSH,Ues,Xt0,Kt,Yt0,Lambda,Eta,Rr,D,E,F,G)
%Calculate Secondary Belt Tenion
%Torsional Spring
%Linear Rate
Fxts = Kt(D).*(Xt0+X2); %Linear spring force (lbf) [Vector]
%Torsional Rate
Fyts = (Lambda(F).*(Yt0(E)+Y2));
%Torsional spring force (lbf*in) [Vector]
%Torque Feedback
FHelix=((Ts.*12)+Fyts)./(2*Rr.*tand(Eta(G,:)));
%Sheave Force
Fs=(FHelix+Fxts)./Beta; %Distributed Clamping Force
Fsx=(2.*Fs.*sin(Beta./2)).*2*tand(Phi/2);
%Total force due to clamping
%Centrifugal Force
Ws=Ws.*(pi/30); %Converting to Radians Per Second
Cs=(1/12).*MBelt .*(Rs.^2).*(Ws.^2);
%Distributed Centrifugal Force
Csx=2.*Cs.*sin(Beta./2); %Total Centrifugal force
%Belt Slack Tenion
Ts0x=(Fsx+Csx)./(1+exp(Ues.*Beta));
%Slack Tension in X Direction
Ts0=zeros(1,n); %Accounting for when belt is greater than or less than 180 degrees
for i=1:1:n
    if Beta(i)<=pi
        Ts0(i)=Ts0x(i)./cos(.5*(pi-Beta(i)));
    elseif Beta(i)>pi
        Ts0(i)=Ts0x(i)./cos(.5*(Beta(i)-pi));
    end
end
%Taught Side Tenion
Ts1=Ts0.*exp(Ues.*Beta);
%Adjusting For Torque Load from Vehicle Weight
Ts0=Ts0-(Tv.*6)./Rs;
Ts1=Ts1+(Tv.*6)./Rs;
%Maximum Torque Tranferable (Without Slip)
%MSMax=(Fsx./cosd(Phi/2)).*Ues.*(Rs./12);
MSMax=(Ts1-Ts0).*(Rs./12); %Maximum torque through shift without slipping (Ft.*lbs)
%Secondary Efficiency
Effs=((MSMax-Ts)./MSMax).*100;
%Secondary Efficiency into gear box
end