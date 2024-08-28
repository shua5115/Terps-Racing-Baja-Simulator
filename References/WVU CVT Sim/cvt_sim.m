%WVU BAJA CVT Analysis
%Sean Skinnr
%Fall 18 - Spring 20
clc
clear
close all
%% Defined Variables
%SCRIPT
n=100; %Number of points during each phase
%Vehicle Mass
Wv=350; %lbs
%Hill Angle
Gamma=20; %Degrees
%Engine
Pee = 7.5; %Engine power engagment (hp)
Pes = 8.5; %Engine power shift (hp)
Wde = 2800; %Desired engagement angular velocity of engine (RPM)
Wds = 3600; %Desired shift Angular Velocity of engine (RPM)
%Tires
RTire = 11; %Tire Radius (in.)
%Gear Box
Rgb=1/7; %Gear Box Reduction Ratio
%BELT
cc = 10; %Center to center (in.)
U = 0.13; %Shigleys Pg 915. Values from Gates Rubber Company
MBelt=0.01579/32.174; %Mass per unit length (Sluggs/in)
BeltH=0.55; %Cross cl Height (in.)
BeltWTop=0.85; %Cross Sectional Width at top (in.)
BeltWBottomn=0.65; %Cross Sectional Width at Bottomn (in.)
Ks = 1; %Service factor was 1.3 from using 1 for testing(Table 17-15 Shigley's)
Kc = 0.965; %Belt parameter (Table 17-16 Shigley's)
Kb = 576; %Belt parameter (Table 17-16 Shigley's)
Nd = 1; %Design factor
%PRIMARY
%Sheaves
Theta = 23; %Belt groove primary (Degrees)
X1max = 0.8; %Total sheave displacement (in.)
Rpmin = 0.75+(BeltH/2); %Minimum pitch radius (in.)
Rpmax = 3-(BeltH/2); %Maximum pitch radius (in.)
%Pressure Spring
Kp=xlsread('CVT Constants 2',1,'N2:Q2');
%Linar spring rates (lbf/in.) [Set: A]
Xp0 = 0.975; %Intalled pretenion (in.)
%Flyweights
Ufly=0.002; %Coefficient of friction between flyweight ramps and flyweight roller bearings
Deltastart=xlsread('CVT Constants 2',1,'D2:F2')'; %Ramp angles (Degrees) [Set: B]
Deltaend=xlsread('CVT Constants 2',1,'D3:F3')';
FlatAngle=xlsread('CVT Constants 2',1,'D9:F9')';
I0=.05; %Flat angle constant. Affects initial shift out point
Mfly=xlsread('CVT Constants 2',1,'D4:I4').*(4/14593.903);
%Flyweight Mass for four arm (Slugs) [Set: C]
Rfly0 = 1.875; %Initial radius (in.) (1.34)
RflyMax=2.345; %Full Shift Radius (in.) (1.815)
Mlink=21.*(4/14593.903); %Link mass for four arms (Sluggs)
%SECONDARY
%Sheaves
Phi = 23; %Belt groove secondary (Degrees)
X2max = 0.8; %Total sheave displacement (in.)
Rsmin = Rpmax*0.9; %Minimum pitch radius (in.)
Rsmax = Rpmin*3.9; %Maximum pitch radius (in.)
%Torsional Spring
Xt0 = 1.44; %Intalled linar Length (in.)
Kt=xlsread('CVT Constants 2',1,'N3:Q3');
%Linar spring rates (lbf/in.) [Set: D]
Yt0=xlsread('CVT Constants 2',1,'D7:I7');
%Intalled twists (Degrees) [Set: E]
Ymax=40; %Total Secondary Twist (Degrees)
Lambda=xlsread('CVT Constants 2',1,'N4:O4');
%Torsional spring rates (in*lbf/Degree.) [Set: F]
%Helix
Etastart = xlsread('CVT Constants 2',1,'D5:H5')'; %Ramp angles (Degrees) [Set: G]
Etaend=xlsread('CVT Constants 2',1,'D6:H6')';
Rr = 1.455; %Ramp radius (in.)
%% Calculated Variables:
%Engine
%Design Hp 66
Hde = Pee*Ks*Nd; %Design power engagement (hp) (Eq 17-19 Shigley's)
Hds = Pes*Ks*Nd; %Design power during shift (hp)
%Vehicle Weight Torque Load
Tv=(Wv*sind(Gamma)*RTire*Rgb)/12;
%ft.*lbs
%Belt
BeltSH=((((BeltWTop-BeltWBottomn)/2)^2)+(BeltH^2))^(.5);
%Side Face Height (in.)
%Flyweights
Rfly=linspace(Rfly0, RflyMax, n);
%Flyweight ramp radius through shift (in.) [Vector]
Rlink=1.625+(0.625.*acos((Rr-1.625)./1.25));
%Flyweight Ramp Angles
Delta=FlatAngle.*ones(length(Deltastart),n);
for i=1:size(Delta,1)
    Delta(i,((I0*n)+1):n)=linspace(Deltastart(i),Deltaend(i),(n-(I0*n)));
    %i=i+1;
end
%Sheaves
Rp = linspace(Rpmin, Rpmax, n); %Primary radius through shift (in.) [Vector]
Rs = linspace(Rsmax, Rsmin, n); %Secondary radius through shift (in.) [Vector]
%Shift Ratio
R=Rs./Rp;
%Shift Distances
X1 = linspace(0, X1max, n); %Primary sheave displacement through shift (in.) [Vector]
X2 = linspace(0, X2max, n); %Secondary sheave displacement through shfit (in.) [Vector]
Y2 = linspace(0, Ymax, n); %Secondary sheave twist through shift (Degrees) [Vector]
%Helix
Eta=Etastart.*ones(length(Etastart),n);
for i=1:size(Eta,1)
    Eta(i,:)=linspace(Eta(i,1),Etaend(i,1),n);
    %i=i+1;
end
%Effective Coefficient of Friction
Uep = U/(sind((Theta)/2)); %Primary effective coefficient of friction
Ues = U/(sind((Phi)/2)); %Secondary effective coefficient of friciton
%% Calculation
%Shift
[TpShift,TsShift,AlphaShift,BetaShift,Tp1Shift,Tp0Shift,CpxShift,...
    FpxShift,FpsShift,FflyxShift,FflyShift,wpcShift,MPMaxShift,FxtsShift,...
    FytsShift,FHelixShift,FsShift,FsxShift,CsxShift,Ts0Shift,Ts1Shift,...
    MSMaxShift,EffsShift,SetupShift] = ...
ShiftV8(n,Hds,Wds,Tv,Rp,Rs,R,X1,X2,Y2,U,Uep,Ues,MBelt,BeltSH,Ufly,cc,Kc,...
    Theta,Phi,Kp,Xp0,Delta,Mfly,Rfly,Mlink,Rlink,Xt0,Kt,Yt0,Lambda,Eta,Rr);
%Engagement Point
[TpEngage,TsEngage,Alpha,Beta,Tp1Engage,Tp0Engage,CpxEngage,FpxEngage,...
    FpsEngage,FflyxEngage,FflyEngage,wpcEngage,MPMaxEngage,FxtsEngage,...
    FytsEngage,FHelixEngage,FsxEngage,CsxEngage,Ts0Engage,Ts1Engage,...
    MSMaxEngage,EffsEngage,SetupEngage] = ...
EngagementV2(n,SetupShift,Hde,Wde,Tv,Rp,Rs,R,U,Uep,Ues,MBelt,Ufly,cc,...
    Theta,Phi,Kp,Xp0,Delta,Mfly,Rfly,Mlink,Rlink,Xt0,Kt,Yt0,Lambda,Eta,Rr);
%Truncating Unacceptable Setups
%Primary
Tp1Shift=Tp1Shift(any(SetupEngage,2),:);
Tp0Shift=Tp0Shift(any(SetupEngage,2),:);
FpxShift=FpxShift(any(SetupEngage,2),:);
FpsShift=FpsShift(any(SetupEngage,2),:);
FflyxShift=FflyxShift(any(SetupEngage,2),:);
FflyShift=FflyShift(any(SetupEngage,2),:);
wpcShift=wpcShift(any(SetupEngage,2),:);
MPMaxShift=MPMaxShift(any(SetupEngage,2),:);
%Secondary
FxtsShift=FxtsShift(any(SetupEngage,2),:);
FytsShift=FytsShift(any(SetupEngage,2),:);
FHelixShift=FHelixShift(any(SetupEngage,2),:);
FsShift=FsShift(any(SetupEngage,2),:);
FsxShift=FsxShift(any(SetupEngage,2),:);
Ts0Shift=Ts0Shift(any(SetupEngage,2),:);
Ts1Shift=Ts1Shift(any(SetupEngage,2),:);
MSMaxShift=MSMaxShift(any(SetupEngage,2),:);
EffsShift=EffsShift(any(SetupEngage,2),:);
%System
SetupEngage=SetupEngage(any(SetupEngage,2),:);
Setup=SetupEngage;
%Vehicle Speed
wpcLow=zeros(size(Setup,1),n);
VLow=zeros(size(Setup,1),n);
VShift=zeros(size(Setup,1),n);
for i=1:size(Setup,1)
    wpcLow(i,:)=linspace(wpcEngage(i),wpcShift(i,1),n);
    VLow(i,:)=(wpcLow(i,:).*Rgb.*60.*2.*pi.*RTire)./(R(1)*63360);
    VShift(i,:)=(wpcShift(i,:).*Rgb.*60.*2.*pi.*RTire)./(R.*63360);
end
% current problem: wpc components have zero size in rows
wpc=[wpcLow wpcShift];
V=[VLow VShift];
Rplot=R(1).*ones(1,n);
Rplot=[Rplot R];
%% Plotting
%% Engine Speed Vs. Shift Ratio
figure
size(Rplot)
size(wpc)
plot(Rplot', wpc')
grid on
set (gca,'xdir','reverse')
set (gca,'FontSize',12)
hold on
%Plot Labels
title('Engine Speed Vs. Shift Ratio')
xlabel('Shift Ratio')
ylabel('Engine Speed [RPM]')
ylim([0 5000])
%Reference Markers
y1=yline(Wds,'--',{'Desired Shift RPM'},'fontsize',12);
y1.LabelVerticalAlignment='middle';
y1.LabelHorizontalAlignment='center';
y1=yline(Wde,'--',{'Desired Engagement RPM'},'fontsize',12);
y1.LabelVerticalAlignment='middle';
y1.LabelHorizontalAlignment='center';
%Reference Lines
x1=xline(1,'--',{'Overdrive'},'fontsize',12);
x1.LabelVerticalAlignment='middle';
x1.LabelHorizontalAlignment='center';
%% Vehicle Speed Vs. Engine Speed
figure
plot(V',wpc')
grid on
hold on
%Plot Axes
set(gca,'FontSize',12)
title('Engine Speed Vs. Vehicle Speed')
xlabel('Vehicle Speed [MPH]')
ylabel('Engine Speed (RPM)')
ylim([0 5000])
xlim([0 max(V(:,2*n))])
%Reference Markers
y1=yline(Wds,'--',{'Desired Shift RPM'},'fontsize',12);
y1.LabelVerticalAlignment='middle';
y1.LabelHorizontalAlignment='center';
y1=yline(Wde,'--',{'Desired Engagement RPM'},'fontsize',12);
y1.LabelVerticalAlignment='middle';
y1.LabelHorizontalAlignment='center';
%% Vehicle Speed Vs. Shift Ratio
figure
plot(Rplot,V')
grid on
set (gca,'xdir','reverse')
set(gca,'FontSize',12)
%Plot Axes
title('Vehicle Speed Vs. Shift Ratio')
xlabel('Shift Ratio')
ylabel('Vehicle Speed [MPH]')
xlim([Rplot(2*n) 4])
%Reference Lines
x1=xline(1,'--',{'Overdrive'},'fontsize',12);
x1.LabelVerticalAlignment='middle';
x1.LabelHorizontalAlignment='center';
%% Primary
%% Torque Received Vs. Maximum Torque Transferrable
figure
plot(R,MPMaxShift)
grid on
hold on
%Plot Axes
set(gca,'FontSize',12)
set(gca,'xdir','reverse')
title('Tranferrable Torque at Primary Through Shift')
xlabel('Shift Ratio')
ylabel('Tranferrable Torque (Ft.*lbs)')
ylim([0 50])
xlim([R(n) R(1)])
%Reference Markers
y1=yline(TpShift,'--',{'Minimum Requirement'},'fontsize',12); %Minimum Torque
y1.LabelVerticalAlignment='middle';
y1.LabelHorizontalAlignment='center';
%% Tension Through Shift
figure
plot(R,Tp0Shift,'--',R,Tp1Shift)
grid on
hold on
%Plot Axes
set(gca,'FontSize',12)
set(gca,'xdir','reverse')
title('Belt Tension Through Shift at Primary')
xlabel('Shift Ratio')
ylabel('Belt Tension (Lbf)')
ylim([-100 200])
xlim([R(n) R(1)])
%Text Box
str={'Dashed Line: Slack Side Tension','Solid Line: Taught Side Tension'};
dim=[0.7 0.6 0.3 0.3];
t=annotation('textbox',dim,'String',str,'FitBoxToText','on','fontsize',12);
%% Clamping Force Through Shift
figure
plot(R,FpxShift)
grid on
%Plot Axes
set(gca,'FontSize',12)
set(gca,'xdir','reverse')
title('Total Clamping Force Through Shift at Primary')
xlabel('Shift Ratio')
ylabel('Clamping Force (Lbf)')
ylim([0 200])
xlim([R(n) R(1)])
%% Secondary
%% Torque Received Vs. Maximum Torque Tranferrable
figure
plot(R,MSMaxShift)
grid on
hold on
%Plot Axes
set(gca,'FontSize',12)
set(gca,'xdir','reverse')
title('Tranferrable Torque at SecondaryThrough Shift')
xlabel('Shift Ratio')
ylabel('Tranferrable Torque (Ft.*lbs)')
ylim([0 80])
xlim([R(n) R(1)])
%Reference Markers
plot(R,TsShift,'--')
text(2.25,20,'-- = MinimumRequirement','FontSize',12);
%label(R,'MinimumRequirement','slope','right')
set(gca,'FontSize',12)
%% Tension Through Shift
figure
plot(R,Ts0Shift,'--',R,Ts1Shift)
grid on
hold on
%Plot Axes
set(gca,'FontSize',12)
set(gca,'xdir','reverse')
title('Belt Tension Through Shift at Secondary')
xlabel('Shift Ratio')
ylabel('Belt Tension (Lbf)')
ylim([-100 200])
xlim([R(n) R(1)])
%Text Box
str={'Dashed Line: Slack Side Tension','Solid Line: Taught Side Tension'};
dim=[0.7 0.6 0.3 0.3];
t=annotation('textbox',dim,'String',str,'FitBoxToText','on','fontsize',12);
%% Clamping Force Through Shift
figure
plot(R,FsxShift)
grid on
%Plot Axes
set(gca,'FontSize',12)
set(gca,'xdir','reverse')
title('Total Clamping Force Through Shift at Secondary')
xlabel('Shift Ratio')
ylabel('Clamping Force (Lbf)')
ylim([0 100])
xlim([R(n) R(1)])
