function []=MainScript(~)
% This is the main script for runing relevant functions

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

Initialization();

%ElectricalPotential(xi_0,pi)
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
hold on

%[startX,startY,startZ]=ProlateSpheroidal2Cartesian_Coor(linspace(xi_0,xi_0,10),linspace(pi/2,pi,10),linspace(0,0,10));
%streamline(X(101:end,101:end),Z(101:end,101:end),EX(101:end,101:end),EZ(101:end,101:end),startX,startZ);
% down mesh for better visualization of electrical field
xi_e=linspace(xi_0,10*xi_0,MeshNum/10);
theta_e=linspace(0,pi,MeshNum/10);
[Xi_e,Theta_e]=meshgrid(xi_e,theta_e);
Phi_e=zeros(size(Xi_e));
[X_e,Y_e,Z_e]=ProlateSpheroidal2Cartesian_Coor(Xi_e,Theta_e,Phi_e);
% Obtain electrical field
EXiComp=ElectricalFieldXiComp(Xi_e,Theta_e);
EThetaComp=ElectricalFieldThetaComp(Xi_e,Theta_e);
[EX,EZ]=ProlateSpheroidal2Cartesian_Comp(Xi_e,Theta_e,EXiComp,EThetaComp);
% quiver plot of electrical field
quiver(X_e,Z_e,EX,EZ,1.5);

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
deltax=1e-6; % x-interval used for calculating
DVal=(LegendreP(l,m,x+deltax)-LegendreP(l,m,x-deltax))/(2*deltax);
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

function [IntVal]=EPIntegral(theta,n)
% This function is an auxiliary function for calculating coefficients of
% electrical potential
global xi_0
Temp=legendre(n,cos(theta));
IntVal=CoverageFunc(theta).*sqrt(sinh(xi_0).^2+sin(theta).^2).*Temp(1,:).*(-sin(theta));
end