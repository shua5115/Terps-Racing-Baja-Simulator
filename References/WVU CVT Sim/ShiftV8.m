function [TpShift,TsShift,AlphaShift,BetaShift,Tp1Shift,Tp0Shift,CpxShift,...
    FpxShift,FpsShift,FflyxShift,FflyShift,wpcShift,MPMaxShift,FxtsShift,...
    FytsShift,FHelixShift,FsShift,FsxShift,CsxShift,Ts0Shift,Ts1Shift,...
    MSMaxShift,EffsShift,SetupShift] = ...
    ShiftV8(n,Hds,Wds,Tv,Rp,Rs,R,X1,X2,Y2,U,Uep,Ues,MBelt,BeltSH,Ufly,cc,Kc...
    ,Theta,Phi,Kp,Xp0,Delta,Mfly,Rfly,Mlink,Rlink,Xt0,Kt,Yt0,Lambda,Eta,Rr)
%Calculate Setup For Shift
FsShift = [];
%% Shift Condition
%Design Power
Hd = Hds; %Design power (hp)
%Angular Velocity
Wp = Wds; %Primary angular velocity (RPM)
Ws=Wp./Rs; %Secondary angular velocity (RPM)
%Torque at Desired RPM
Tp=(Hd.*5252)./Wp; %Primary Torque (ft.*lbs)
Ts=Tp.*R; %Torque at secondary through shift (ft.*lbs)
%% Belt Calculation
[Alpha,Beta] = BeltV2(n,Rp,Rs,cc);
%% Preallocation
%Primary
szx = length(Kp)*size(Delta,1)*length(Mfly)*length(Kt)*length(Yt0)*length(Lambda)*size(Eta,1);
szy = n;
Tp1Shift=zeros(szx, szy);
Tp0Shift=zeros(szx, szy);
CpxShift=zeros(szx, szy);
FpxShift=zeros(szx, szy);
FpsShift=zeros(szx, szy);
FflyxShift=zeros(szx, szy);
FflyShift=zeros(szx, szy);
wpcShift=zeros(szx, szy);
MPMaxShift=zeros(szx, szy);
%Secondary
FxtsShift=zeros(szx, szy);
FytsShift=zeros(szx, szy);
FHelixShift=zeros(szx, szy);
CsxShift=zeros(szx, szy);
FsxShift=zeros(szx, szy);
Ts0Shift=zeros(szx, szy);
Ts1Shift=zeros(szx, szy);

MSMaxShift=zeros(szx, szy);
EffsShift=zeros(szx, szy);
%System
szxx = length(Kp)*size(Delta,1)*length(Mfly)*length(Kt)*length(Yt0)*length(Lambda)*size(Eta,1);
SetupShift=zeros(szxx,7);
%% Primary and Secondary Calculation
i=0;
for A=1:length(Kp)
    for B=1:size(Delta,1)
        for C=1:length(Mfly)
            for D=1:length(Kt)
                for E=1:length(Yt0)
                    for F=1:length(Lambda)
                        for G=1:size(Eta,1)
[Fxts,Fyts,FHelix,Fs,Fsx,Csx,Ts0,Ts1,MSMax,Effs] =...
    SecondaryV3(n,MBelt,Ws,Tv,Ts,X2,Y2,Rs,Phi,Beta,BeltSH,Ues,Xt0,Kt,Yt0,Lambda,Eta,Rr,D,E,F,G);
[Tp1,Tp0,Cpx,Fpx,Fps,Fflyx,Ffly,wpc,MPMax] =...
    PrimaryV3(n,Tp,MBelt,Ts1,Wp,X1,Rp,Alpha,BeltSH,Uep,Ufly,Theta,Kp,Xp0,Delta,Mfly,Rfly,Mlink,Rlink,A,B,C);
                            if any(wpc>3400 & wpc<3800 & MPMax>=0.8*Tp & MSMax>=0.8*Ts, 'all') %& MPMax<=1.3.*Tp
                                i=i+1;
                                %Compiling Data Points
                                %Primary
                                Tp1Shift(i,:)=Tp1;
                                Tp0Shift(i,:)=Tp0;
                                CpxShift(i,:)=Cpx;
                                FpxShift(i,:)=Fpx;
                                FpsShift(i,:)=Fps;
                                FflyxShift(i,:)=Fflyx;
                                FflyShift(i,:)=Ffly;
                                wpcShift(i,:)=wpc;
                                MPMaxShift(i,:)=MPMax;
                                %Secondary
                                FxtsShift(i,:)=Fxts;
                                FytsShift(i,:)=Fyts;
                                FHelixShift(i,:)=FHelix;
                                FsShift(i,:)=Fs;
                                FsxShift(i,:)=Fsx;
                                CsxShift(i,:)=Csx;
                                Ts0Shift(i,:)=Ts0;
                                Ts1Shift(i,:)=Ts1;
                                MSMaxShift(i,:)=MSMax;
                                EffsShift(i,:)=Effs;
                                %System
                                SetupShift(i,1)=A;
                                SetupShift(i,2)=B;
                                SetupShift(i,3)=C;
                                SetupShift(i,4)=D;
                                SetupShift(i,5)=E;
                                SetupShift(i,6)=F;
                                SetupShift(i,7)=G;
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Truncate unacceptable combination
%Primary
Tp1Shift=Tp1Shift(any(Tp1Shift,2),:);
Tp0Shift=Tp0Shift(any(Tp0Shift,2),:);
CpxShift=CpxShift(any(FpxShift,2),:);
FpxShift=FpxShift(any(FpxShift,2),:);
FpsShift=FpsShift(any(FpsShift,2),:);
FflyxShift=FflyxShift(any(FflyxShift,2),:);
FflyShift=FflyShift(any(FflyShift,2),:);
wpcShift=wpcShift(any(wpcShift,2),:);
MPMaxShift=MPMaxShift(any(MPMaxShift,2),:);
%Secondary
FxtsShift=FxtsShift(any(FxtsShift,2),:);
FytsShift=FytsShift(any(FytsShift,2),:);
FHelixShift=FHelixShift(any(FHelixShift,2),:);
FsShift=FsShift(any(FsShift,2),:);
FsxShift=FsxShift(any(FsxShift,2),:);
CsxShift=CsxShift(any(CsxShift,2),:);
Ts0Shift=Ts0Shift(any(Ts0Shift,2),:);
Ts1Shift=Ts1Shift(any(Ts1Shift,2),:);
MSMaxShift=MSMaxShift(any(MSMaxShift,2),:);
EffsShift=EffsShift(any(EffsShift,2),:);
%System
SetupShift=SetupShift(any(SetupShift,2),:); % problem: this filters out everything
%% Renaming
TpShift=Tp;
TsShift=Ts;
AlphaShift=Alpha;
BetaShift=Beta;
end