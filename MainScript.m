function []=MainScript(~)
% This is the main script for runing relevant functions

% Set the value of global parameters
ParameterSetting();
% Initialization
Initialization();

DrawElectricalPotential()
%hold on
%DrawElectricalField()
%hold off

figure
%DrawStreamFunction()
%DrawVelocityFieldMag()
hold on
%DrawVelocityField()
i=30;
theta=linspace(0,pi,1000);
%plot(x,LegendreP(i,0,x))
%plot(x,DLegendreP(i,0,x))
%plot(x,DDLegendreP(i,0,x))
%legend('0','1','2')
xi=zeros(size(theta));
%plot(theta,VelocityXiCompSeries(xi,theta));
%plot(theta,VelocityThetaCompSeries(xi,theta));
%global EPTruncNum
%global EPCoefList
%plot(0:1:EPTruncNum,EPCoefList)
end

function []=Initialization()
% This part is used to initialize useful data
global EPTruncNum
global EPCoefList

% Prepare EP coefficients
for i=0:1:EPTruncNum
    EPCoefList(i+1)=CoefElectricalPotential(i);
end

% Prepare LegendreP

% Prepare DLegendreP

% Prepare LegendreQ

% Prepare DLegendreQ
end

function []=ParameterSetting()
%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
global a% Length of Semi-Major Axis Unit:um
a=6;
global b % Length of Semi-Minor Axis Unit:um
b=1;
global c; % Length of Focal Length
c=sqrt(a^2-b^2);
global xi_0 % Parameter \xi_0 which defines motor surface at prolate spheroidal coordinate
xi_0=acosh(a/c);
global epsilon_0 % Electrical Permivittivity
epsilon_0=1;
global sigma_pos % Surface Charge Density at Positive Side
sigma_pos=1;
global sigma_neg % Surface Charge Density at Negative Side
sigma_neg=-1;
global BiPosTheta % Theta which characterizes the boundary of positive side and negative side
BiPosTheta=pi/2;
global Sharpness % Control the width of transition zone between positive side and negative side
Sharpness=1e-5;
global LegendreQTruncNum % Control the number of terms which used for calculating the second kind of Legendre function
LegendreQTruncNum=10;
global EPTruncNum % Control the number of terms which used for calculating electrical potential
EPTruncNum=30;
global EPCoefList % Store coefficients appeared in electrical potential for repeating usage
EPCoefList=zeros(1,EPTruncNum+1);
global MeshNum
MeshNum=200;
global U_motor % Define the velocity of motor speed (\mu m /s)
U_motor=0;
global mu_motor % Define the mobility of motor (unit)
mu_motor=1;
end

function []=DrawElectricalPotential()
global xi_0
global MeshNum

xi=linspace(xi_0,10*xi_0,MeshNum);
theta=linspace(0,pi,MeshNum);
[Xi,Theta]=meshgrid(xi,theta);
% Because of the axis symmetry, fixed Phi=0
Phi=zeros(size(Xi));
% Convert to Cartesian coordinates
[X,Y,Z]=ProlateSpheroidal2Cartesian_Coor(Xi,Theta,Phi);
EPVal=ElectricalPotential(Xi,Theta);
% contour plot of electrical potential
contourf(X,Z,EPVal)
end

function []=DrawElectricalField()
global xi_0
global MeshNum
%[startX,startY,startZ]=ProlateSpheroidal2Cartesian_Coor(linspace(xi_0,xi_0,10),linspace(pi/2,pi,10),linspace(0,0,10));
%streamline(X(101:end,101:end),Z(101:end,101:end),EX(101:end,101:end),EZ(101:end,101:end),startX,startZ);
% down mesh for better visualization of electrical field
xi=linspace(xi_0,10*xi_0,MeshNum/10);
theta=linspace(0,pi,MeshNum/10);
[Xi,Theta]=meshgrid(xi,theta);
Phi=zeros(size(Xi));
[X,Y,Z]=ProlateSpheroidal2Cartesian_Coor(Xi,Theta,Phi);
% Obtain electrical field
EXiComp=ElectricalFieldXiComp(Xi,Theta);
EThetaComp=ElectricalFieldThetaComp(Xi,Theta);
[EX,EZ]=ProlateSpheroidal2Cartesian_Comp(Xi,Theta,EXiComp,EThetaComp);
% quiver plot of electrical field
quiver(X,Z,EX,EZ,1.5);
end

function []=DrawStreamFunction()
global xi_0
global MeshNum

xi=linspace(xi_0,10*xi_0,MeshNum*3);
theta=linspace(0,pi,MeshNum*2);
[Xi,Theta]=meshgrid(xi,theta);
% Because of the axis symmetry, fixed Phi=0
Phi=zeros(size(Xi));
% Convert to Cartesian coordinates
[X,Y,Z]=ProlateSpheroidal2Cartesian_Coor(Xi,Theta,Phi);
EThetaComp=ElectricalFieldThetaComp(Xi,Theta);
StreamFuncVal=StokesStreamFunc(Xi,Theta,EThetaComp);
% contour plot of electrical potential
contourf(X,Z,StreamFuncVal)
end

function []=DrawVelocityFieldMag()
global xi_0
global MeshNum

xi=linspace(xi_0,10*xi_0,MeshNum*3);
theta=linspace(0,pi,MeshNum*2);
[Xi,Theta]=meshgrid(xi,theta);
% Because of the axis symmetry, fixed Phi=0
Phi=zeros(size(Xi));
% Convert to Cartesian coordinates
[X,Y,Z]=ProlateSpheroidal2Cartesian_Coor(Xi,Theta,Phi);
% Obtain velocity field
VXiComp=VelocityFieldXiComp(Xi,Theta);
VThetaComp=VelocityFieldThetaComp(Xi,Theta);
VelocityMag=sqrt(VXiComp.^2+VThetaComp.^2);
% contour plot of electrical potential
VelocityMag=ExcludeOutliers(VelocityMag,93,0.001);
contourf(X,Z,VelocityMag);
end

function []=DrawVelocityField()
global xi_0
global MeshNum
xi=linspace(xi_0,10*xi_0,MeshNum/2);
theta=linspace(0+0.01,pi-0.01,MeshNum/2);
[Xi,Theta]=meshgrid(xi,theta);
Phi=zeros(size(Xi));
[X,Y,Z]=ProlateSpheroidal2Cartesian_Coor(Xi,Theta,Phi);
% Obtain velocity field
VXiComp=VelocityFieldXiComp(Xi,Theta);
VThetaComp=VelocityFieldThetaComp(Xi,Theta);
[VX,VZ]=ProlateSpheroidal2Cartesian_Comp(Xi,Theta,VXiComp,VThetaComp);
% clear infinity value
% quiver plot of electrical field
VX=ExcludeOutliers(VX,95,0.01);
VZ=ExcludeOutliers(VZ,95,0.01);
VX_norm=VX./sqrt(VX.^2+VZ.^2);
VZ_norm=VZ./sqrt(VX.^2+VZ.^2);
quiver(X,Z,VX_norm,VZ_norm,1);
%[startX,startY,startZ]=ProlateSpheroidal2Cartesian_Coor(linspace(xi_0,xi_0,10),linspace(pi/2,pi,10),linspace(0,0,10));
%streamline(X,Z,VX_norm,VZ_norm,startX,startZ);
end

function [X,Y,Z]=ProlateSpheroidal2Cartesian_Coor(Xi,Theta,Phi)
% Transform coordinates from prolate spheroidal coordinates system to
% cartesian coordinates system
global c
X=c*sinh(Xi).*sin(Theta).*cos(Phi);
Y=c*sinh(Xi).*sin(Theta).*sin(Phi);
Z=c*cosh(Xi).*cos(Theta);
end

function [XComp,ZComp]=ProlateSpheroidal2Cartesian_Comp(Xi,Theta,XiComp,ThetaComp)
% Transform Xi component and Theta Component in prolatespheroidal
% coordinates system to cartesian coordinates system
global c
% angle which characterizes the relation between two coordinate
% corresponding x coordinates value
%X=c.*sinh(Xi).*sin(Theta);
% prepare intermediate value
%TVal=(cosh(Xi)./(c*sinh(Xi).^2)).*(X./sqrt(1-(X.^2)./(c^2*sinh(Xi).^2)));
TVal=(cosh(Xi)./sinh(Xi)).*(sin(Theta)./abs(cos(Theta)));
% do inverse arctan operation, specially, here we want it returned value
% range from 0 to pi
Angle=atan(TVal);
Angle(Theta>pi/2)=pi-Angle(Theta>pi/2);

% Then, transform XiComp and ThetaComp into XComp and ZComp according to
% the value of Angle
XComp=XiComp.*sin(Angle)+ThetaComp.*cos(Angle);
ZComp=XiComp.*cos(Angle)-ThetaComp.*sin(Angle);
end

function [EVal]=ElectricalFieldXiComp(xi,theta)
global c
global EPCoefList
global EPTruncNum
global LegendreQTruncNum
EVal=0;

for i=0:1:EPTruncNum
% Loop through each term
        EVal=EVal+sinh(xi)./(c*sqrt(sinh(xi).^2+sin(theta).^2))*EPCoefList(i+1).*DLegendreQ(i,cosh(xi),LegendreQTruncNum).*LegendreP(i,0,cos(theta));
end

end

function [EVal]=ElectricalFieldThetaComp(xi,theta)
% Used for calculating Theta component of electrical field
global c
global EPCoefList
global EPTruncNum
global LegendreQTruncNum

EVal=0;

for i=0:1:EPTruncNum
    % loop till truncted term
        EVal=EVal-sin(theta)./(c*sqrt(sinh(xi).^2+sin(theta).^2))*EPCoefList(i+1).*LegendreQ(i,cosh(xi),LegendreQTruncNum).*DLegendreP(i,0,cos(theta));
end
end

function [EPVal]=ElectricalPotential(xi,theta)
% This function is used for calculating the electrical potential
global EPTruncNum
global LegendreQTruncNum
global EPCoefList

EPVal=0;
for i=0:1:EPTruncNum
        EPVal=EPVal+EPCoefList(i+1)*LegendreQ(i,cosh(xi),LegendreQTruncNum).*LegendreP(i,0,cos(theta));
end
end

function [ ChargeDensity ] = CoverageFunc( theta)
%COVERAGEFUNC Coverage function used to describe the charge distribution at
%motor surface
global sigma_pos
global sigma_neg
global BiPosTheta
global Sharpness
ChargeDensity=sigma_pos*tanh((theta-BiPosTheta)/Sharpness);
end

function [Val] = LegendreQ(n,x,TruncNum)
PreFactor=1;
Factor=1;
Val=0;
for i=1:1:TruncNum
% Loop till the truncation limit
Val=Val+Factor*x.^(-(n+2*i-1));
Factor=Factor*(n+i)*(n+i+1)/((2*i)*(2*n+2*i+1));
end
for i=1:1:n
    PreFactor=PreFactor*i/(2*i-1);
end
Val=Val*PreFactor;
end

function [DVal] = DLegendreQ(n,x,TruncNum)
PreFactor=1;
Factor=1;
DVal=0;
for i=1:1:TruncNum
% Loop till the truncation limit
DVal=DVal+Factor*(-(n+2*i-1))*x.^(-(n+2*i));
Factor=Factor*(n+i)*(n+i+1)/((2*i)*(2*n+2*i+1));
end
for i=1:1:n
    PreFactor=PreFactor*i/(2*i-1);
end
DVal=DVal*PreFactor;
end

function [DVal] = DLegendreP(l,m,x)
% This function is used for calculating derivatives of legendre polynomials
% simple central difference method
deltax=1e-5; % x-interval used for calculating
DVal=(LegendreP(l,m,x+deltax)-LegendreP(l,m,x-deltax))/(2*deltax);
end

function [DDVal] = DDLegendreP(l,m,x)
% This function is used for calculating second derivatives of legendre
% polynoials
deltax=1e-2;
DDVal=(DLegendreP(l,m,x+deltax)-DLegendreP(l,m,x-deltax))/(2*deltax);
end

function [ Coef ] = CoefElectricalPotential( n )
%COEFELECTRICALPOTENTIAL Summary of this function goes here
%   Detailed explanation goes here
global xi_0
global c
global epsilon_0
global LegendreQTruncNum
Coef=(2*n+1)/2*(1/DLegendreQ(n,cosh(xi_0),LegendreQTruncNum))*(c/(epsilon_0*sinh(xi_0)))*integral(@(theta)EPIntegral(theta,n),pi,0);
end

function [Val] = StokesStreamFunc(Xi,Theta,ETheta)
% This function is used for calculating the value of stokes stream function
% global parameters
global c
global xi_0
global U_motor
global mu_motor
% auxiliary variables
lambda_0=cosh(xi_0);
Lambda=cosh(Xi);
Kappa=cos(Theta);
Omega=c*sqrt(Lambda.^2-1).*sqrt(1-Kappa.^2);
sigma_1=-lambda_0+lambda_0^2*acoth(lambda_0)+acoth(lambda_0);
sigma_2=-lambda_0+lambda_0^2*acoth(lambda_0)-acoth(lambda_0);
% calculation
% Hydrodynamic part
HydroUpper=Lambda./(Lambda.^2-1)-(lambda_0^2+1)/(lambda_0^2-1)*acoth(Lambda);
HydroLower=lambda_0/(lambda_0^2-1)-(lambda_0^2+1)/(lambda_0^2-1)*acoth(lambda_0);
HydrodynamicPart=-0.5*U_motor*Omega.^2.*(HydroUpper./HydroLower);
% Electrostatic part
ElectroFactor=sqrt(lambda_0^2-Kappa.^2)./((lambda_0^2-1)*sqrt(1-Kappa.^2));
ElectroPart=-0.5*sigma_2/sigma_1*ElectroFactor.*(Lambda-2*lambda_0)*mu_motor.*ETheta;
% Add these two part together
Val=HydrodynamicPart+ElectroPart;
end

function [VXiVal] = VelocityFieldXiComp(Xi,Theta)
% This function is used for calculating the xi component of velocity field
% global parameters
global c
global xi_0
global U_motor
global mu_motor
% auxiliary variables
Lambda=cosh(Xi);
lambda_0=cosh(xi_0);
Kappa=cos(Theta);
sigma_1=-lambda_0+lambda_0^2*acoth(lambda_0)+acoth(lambda_0);
sigma_2=-lambda_0+lambda_0^2*acoth(lambda_0)-acoth(lambda_0);
% Calculation
% Hydrodynamic Part
HydroFactor=U_motor*Kappa.*sqrt(Lambda.^2-1)./sqrt(Lambda.^2-Kappa.^2);
HydroUpper=Lambda./(Lambda.^2-1)-(lambda_0.^2+1)/(lambda_0.^2-1)*acoth(Lambda);
HydroLower=lambda_0/(lambda_0^2-1)-(lambda_0.^2+1)/(lambda_0.^2-1)*acoth(lambda_0);
HydrodynamicPart=HydroFactor.*HydroUpper./HydroLower;
% Electrostatic Part
ElectroFactor1=0.5*(Lambda-2*lambda_0)./(c^2*sqrt(Lambda.^2-1).*sqrt(Lambda.^2-Kappa.^2));
ElectroFactor2=1/(lambda_0^2-1)*(sigma_2/sigma_1)*(mu_motor/c);
ElectroSeriesVal=VelocityXiCompSeries(Xi,Theta);
ElectroPart=ElectroFactor1.*ElectroFactor2.*ElectroSeriesVal;
% Combine together
VXiVal=HydrodynamicPart+ElectroPart;
end

function [VThetaVal] = VelocityFieldThetaComp(Xi,Theta)
% This function is used for calculating the xi component of velocity field
% global parameters
global c
global xi_0
global U_motor
global mu_motor
% auxiliary variables
Lambda=cosh(Xi);
lambda_0=cosh(xi_0);
Kappa=cos(Theta);
sigma_1=-lambda_0+lambda_0^2*acoth(lambda_0)+acoth(lambda_0);
sigma_2=-lambda_0+lambda_0^2*acoth(lambda_0)-acoth(lambda_0);
% Calculation
% Hydrodynamic Part
HydroFactor=-U_motor*sqrt(1-Kappa.^2)./(2*sqrt(Lambda.^2-Kappa.^2));
HydroUpper=1-((lambda_0^2+1)/(lambda_0^2-1)*(2*Lambda.*acoth(Lambda)-1));
HydroLower=lambda_0/(lambda_0^2-1)-(lambda_0^2+1)/(lambda_0^2-1)*acoth(lambda_0);
HydrodynamicPart=HydroFactor.*HydroUpper./HydroLower;
% Electrostatic Part
ElectroFactor1=1./(0.5*(c^2.*sqrt(1-Kappa.^2).*sqrt(Lambda.^2-Kappa.^2)));
ElectroFactor2=(1/(lambda_0^2-1))*(sigma_2/sigma_1)*(mu_motor/c);
ElectroSeriesVal=VelocityThetaCompSeries(Xi,Theta);
ElectroPart=ElectroFactor1.*ElectroFactor2.*ElectroSeriesVal;
% Combine together
VThetaVal=HydrodynamicPart+ElectroPart;
end

function [Val] = VelocityThetaCompSeries(Xi,Theta)
global EPTruncNum
global LegendreQTruncNum
global EPCoefList
global xi_0
Val=0;
for i=0:1:EPTruncNum
    Val=Val+EPCoefList(i+1)*LegendreQ(i,cosh(xi_0),LegendreQTruncNum)*DLegendreP(i,0,cos(Theta));
end
end

function [Val] = VelocityXiCompSeries(Xi,Theta)
global EPTruncNum
global LegendreQTruncNum
global EPCoefList
global xi_0
Val=0;
for i=0:1:EPTruncNum
    Val=Val+EPCoefList(i+1)*LegendreQ(i,cosh(xi_0),LegendreQTruncNum).*DDLegendreP(i,0,cos(Theta));
end
end

function [IntVal]=EPIntegral(theta,n)
% This function is an auxiliary function for calculating coefficients of
% electrical potential
global xi_0
Temp=legendre(n,cos(theta));
IntVal=CoverageFunc(theta).*sqrt(sinh(xi_0).^2+sin(theta).^2).*Temp(1,:).*(-sin(theta));
end

function [CleanedData] = ExcludeOutliers(RawData,percentileUp,percentileDown)
% This function is used for clearing very large and small value
UpperLim=prctile(RawData(:),percentileUp);
LowerLim=prctile(RawData(:),percentileDown);
RawData((RawData>UpperLim) | (RawData<LowerLim))=NaN;
CleanedData=RawData;
end