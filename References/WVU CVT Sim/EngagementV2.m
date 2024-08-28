function [TpEngage,TsEngage,Alpha,Beta,Tp1Engage,Tp0Engage,CpxEngage,...
    FpxEngage,FpsEngage,FflyxEngage,FflyEngage,wpcEngage,MPMaxEngage,...
    FxtsEngage,FytsEngage,FHelixEngage,FsxEngage,CsxEngage,Ts0Engage,...
    Ts1Engage,MSMaxEngage,EffsEngage,SetupEngage] = ...
    EngagementV2(n,SetupShift,Hde,Wde,Tv,Rp,Rs,R,U,Uep,Ues,MBelt,Ufly,cc,Theta,...
    Phi,Kp,Xp0,Delta,Mfly,Rfly,Mlink,Rlink,Xt0,Kt,Yt0,Lambda,Eta,Rr)
%Calculate Setup For Optimal Engagement Point
%% Engagement Condition
%Design Power
Hd = Hde; %Design power at engagement (hp)
%Sheaves
Rp=Rp(1); %Primary radius during engagement
Rs=Rs(1); %Secondary radius during engagement
%Angular Velocity
Wp = Wde; %Primary angular velocity (RPM) *Assuming minimal speed differential*
Ws=Wp./Rs; %Secondary angular velocity (RPM)
%Flyweights
Rfly=Rfly(1);
Rlink=Rlink(1);
%Torque at Desired RPM
Tp=(Hd.*5252)./Wp; %Primary Torque (ft.*lbs)
Ts=Tp.*R(1); %Torque at secondary through shift (ft.*lbs)
%% Preallocation
%Primary
Tp1Engage=zeros(size(SetupShift,1),1);
Tp0Engage=zeros(size(SetupShift,1),1);
CpxEngage=zeros(size(SetupShift,1),1);
FpxEngage=zeros(size(SetupShift,1),1);
FpsEngage=zeros(size(SetupShift,1),1);
FflyxEngage=zeros(size(SetupShift,1),1);
FflyEngage=zeros(size(SetupShift,1),1);
wpcEngage=zeros(size(SetupShift,1),1);
MPMaxEngage=zeros(size(SetupShift,1),1);
%Secondary
FxtsEngage=zeros(size(SetupShift,1),1);
FytsEngage=zeros(size(SetupShift,1),1);
FHelixEngage=zeros(size(SetupShift,1),1);
CsxEngage=zeros(size(SetupShift,1),1);
FsxEngage=zeros(size(SetupShift,1),1);
Ts0Engage=zeros(size(SetupShift,1),1);
Ts1Engage=zeros(size(SetupShift,1),1);
MSMaxEngage=zeros(size(SetupShift,1),1);
EffsEngage=zeros(size(SetupShift,1),1);
%System
SetupEngage=zeros(size(SetupShift,1),7);
KpEngage=zeros(1,size(SetupShift,1));
DeltaEngage=zeros(1,size(SetupShift,1));
MflyEngage=zeros(1,size(SetupShift,1));
KtEngage=zeros(1,size(SetupShift,1));
Yt0Engage=zeros(1,size(SetupShift,1));
LambdaEngage=zeros(1,size(SetupShift,1));
EtaEngage=zeros(1,size(SetupShift,1));
%% Possible Engagement Setups
for i=1:size(SetupShift,1)
    KpEngage(i)=Kp(SetupShift(i,1))';
    DeltaEngage(i)=Delta(SetupShift(i,2))';
    MflyEngage(i)=Mfly(SetupShift(i,3))';
    KtEngage(i)=Kt(SetupShift(i,4))';
    Yt0Engage(i)=Yt0(SetupShift(i,5))';

    LambdaEngage(i)=Lambda(SetupShift(i,6))';
    EtaEngage(i)=Eta(SetupShift(i,7))';
end
%% Belt Calculations
Alpha = (pi-(2*asin((Rs-Rp)/cc)));
Beta= (pi+(2*asin((Rs-Rp)/cc)));
%% Primray and Secondary Calculation
for i=1:size(SetupShift,1)
    Kp=KpEngage(i);
    Delta=DeltaEngage(i);
    Mfly=MflyEngage(i);
    Kt=KtEngage(i);
    Yt0=Yt0Engage(i);
    Lambda=LambdaEngage(i);
    Eta=EtaEngage(i);

    %% Secondary
    %Torsional Spring
    %Linear Rate
    Fxts = Kt.*(Xt0); %Linear spring force (lbf) [Vector]

    %Torsional Rate
    Fyts = (Lambda.*(Yt0)); %Torsional spring force (lbf*in) [Vector]

    %Torque Feedback

    FHelix=((Ts.*12)+Fyts)./(2*Rr.*tand(Eta));

    %Sheave Force
    Fs=(FHelix+Fxts)./Beta; %Distributed Clamping Force

    Fsx=(2.*Fs.*sin(Beta./2)).*2*tand(Phi/2);
    %Total force due to clamping

    %Centrifugal Force
    Ws=Ws.*(pi/30); %Converting to Radians Per Second
    Cs=(1/12).*MBelt.*(Rs.^2).*(Ws.^2);
    %Distributed Centrifugal Force
    Csx=2.*Cs.*sin(Beta./2);
    %Total Centrifugal force

    %Belt Slack Tenion
    Ts0x=(Fsx+Csx)./(1+exp(Ues.*Beta));

    if Beta<=pi
        Ts0=Ts0x./cos(.5*(pi-Beta));
    elseif Beta>pi
        Ts0=Ts0x./cos(.5*(Beta-pi));
    end


    %Taught Side Tenion
    Ts1=Ts0.*exp(Ues.*Beta);

    %Adjusting For Torque Load from Vehicle Weight
    Ts0=Ts0-(Tv.*6)./Rs;
    Ts1=Ts1+(Tv.*6)./Rs;

    %Maximum Torque Tranferable (Without Slip)
    MSMax=(Ts1-Ts0).*(Rs./12); %Maximum torque through shift without slipping (Ft.*lbs)

    %Secondary Efficiency
    Effs=(Ts./MSMax).*100; %Secondary Efficiency into gear box

    %% Primary
    %Taught side tenion
    Tp1=Ts1;

    %Slack side tenion
    Tp0=Tp1./exp(Uep.*Alpha);

    %Accounting for Motor Torque
    Tp1=Tp1+((Tp.*6)./Rp);
    Tp0=Tp0-((Tp.*6)./Rp);

    %Tension Component in X Direction

    if Alpha<=pi
        Tp0x=Tp0.*cos(.5*(pi-Alpha));
    elseif Alpha>pi
        Tp0x=Tp0.*cos(.5*(Alpha-pi));
    end

    if Alpha<=pi
        Tp1x=Tp1.*cos(.5*(pi-Alpha));
    elseif Alpha>pi
        Tp1x=Tp1.*cos(.5*(Alpha-pi));
    end


    %Centrifugal Force
    Wp=Wp*(pi/30);
    Cp=(1/12).*MBelt.*(Rp.^2).*(Wp.^2); %Distributed Centrifugal Force
    Cpx=2.*Cp.*sin(Alpha./2); %Total Centrifugal Force

    %Required side force
    Fpx=(-Cpx); %Total clamping force

    Fp=(Fpx./(2.*sin(Alpha./2)))./(2.*tand(Theta/2)); %Distributed clamping force

    %Pressure Spring
    Fps=Kp.*(Xp0); %Pressure spring force vector (lbf)

    %Required Flyweight force
    Fflyx=Fp+Fps;
    Ffly=(Fflyx./tand(Delta))-(Fflyx./tan(acos((Rfly-1.625)./1.25)));

    %Required engine speed
    wpc=(12.*Ffly./((Mfly.*Rfly)+(Mlink.*Rlink))).^.5;
    wpc=(wpc.*30)./(pi); %Rev/Min

    %Maximum Torque Tranferable (WithoutSlip)
    MPMax=(Tp1-Tp0).*(Rp./12); %Maximum torque through shift without slipping (Ft.*lbs)

    if wpc>=2100 %& wpc<=3200 & MSMax>=Ts

        %Compiling Data Points
        %Primary
        Tp1Engage(i,:)=Tp1;
        Tp0Engage(i,:)=Tp0;
        CpxEngage(i,:)=Cpx;
        FpxEngage(i,:)=Fpx;
        FpsEngage(i,:)=Fps;
        FflyxEngage(i,:)=Fflyx;
        FflyEngage(i,:)=Ffly;
        wpcEngage(i,:)=wpc;
        MPMaxEngage(i,:)=MPMax;

        %Secondary
        FxtsEngage(i,:)=Fxts;
        FytsEngage(i,:)=Fyts;
        FHelixEngage(i,:)=FHelix;
        CsxEngage(i,:)=Csx;
        FsxEngage(i,:)=Fsx;
        Ts0Engage(i,:)=Ts0;
        Ts1Engage(i,:)=Ts1;
        MSMaxEngage(i,:)=MSMax;
        EffsEngage(i,:)=Effs;

        %System
        SetupEngage(i,1)=SetupShift(i,1);
        SetupEngage(i,2)=SetupShift(i,2);
        SetupEngage(i,3)=SetupShift(i,3);
        SetupEngage(i,4)=SetupShift(i,4);
        SetupEngage(i,5)=SetupShift(i,5);
        SetupEngage(i,6)=SetupShift(i,6);
        SetupEngage(i,7)=SetupShift(i,7);

        i=i+1;
    end
end
%% Truncate unacceptable combination
%Primary
Tp1Engage=Tp1Engage(any(Tp1Engage,2),:);
Tp0Engage=Tp0Engage(any(Tp0Engage,2),:);
CpxEngage=CpxEngage(any(CpxEngage,2),:);
FpxEngage=FpxEngage(any(FpxEngage,2),:);
FpsEngage=FpsEngage(any(FpsEngage,2),:);
FflyxEngage=FflyxEngage(any(FflyxEngage,2),:);
FflyEngage=FflyEngage(any(FflyEngage,2),:);
wpcEngage=wpcEngage(any(wpcEngage,2),:);
MPMaxEngage=MPMaxEngage(any(MPMaxEngage,2),:);
%Secondary
FxtsEngage=FxtsEngage(any(FxtsEngage,2),:);
FytsEngage=FytsEngage(any(FytsEngage,2),:);
FHelixEngage=FHelixEngage(any(FHelixEngage,2),:);
CsxEngage=CsxEngage(any(CsxEngage,2),:);
FsxEngage=FsxEngage(any(FsxEngage,2),:);
Ts0Engage=Ts0Engage(any(Ts0Engage,2),:);
Ts1Engage=Ts1Engage(any(Ts1Engage,2),:);
MSMaxEngage=MSMaxEngage(any(MSMaxEngage,2),:);
EffsEngage=EffsEngage(any(EffsEngage,2),:);
%% Renaming
TpEngage=Tp;
TsEngage=Ts;
AlphaEngage=Alpha;
BetaEngage=Beta;
end
