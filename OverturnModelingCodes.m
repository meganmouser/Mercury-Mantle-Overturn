%Modeling codes for: "On the potential for cumulate mantle overturn in
%Mercury" by M. D. Mouser and N. Dygert

%Instructions: to perform analysis, input values into empty brackets
%corresponding to specific parameters. See main text for demonstration of
%analysis and for more information on parameter selection. 
%%
%Density model (Section 2)
clear all
clc
%determining the pressure and temperature relationship on the molar volume 
V0=[];%inital molar volume
Kprimet=[];%pressure derivitative
Kt298=[];%reference Kt value
dKtdT=[];%P-T relationship
temp=[];%in kelvin
Kt=Kt298+(dKtdT.*(temp-298));%bulk modulus scaling relationship (eq. 2)
P=[];%pressure in GPa

%equation of state (eq. 1)
f=@(Vf) (3/2).*Kt.*(((V0/Vf).^(7/3))-((V0/Vf).^(5/3))).*(1-(3/4).*(4-Kprimet).*(((V0/Vf).^(2/3))-1))-P

fzero(f,300)

%density (d, molar mass-molar volume realtionship)
Vf=[];%input Vf value from fzero
mm=[];%molar mass of material
d=mm./Vf; %density in g/cm3

%%
%Heat production calculation (Section 4)
clear all 
clc
%present day concentrations for isotopes x, y, z (per kg)
x=[];%isotope x
y=[];%isotope y
z=[];%isotope z
%half life (seconds)
hlx=[];%half life for isotope x
hly=[];%half life for isotope y
hlz=[];%half life for isotope z
%decay constant (lambda)
lx=[];%for isotope x
ly=[];%for isotope y
lz=[];%for isotope z
%4.5 Ga in seconds
tGa=1.4193e17;
%number of half lives since decay started
nx=tGa./hlx;%isotope x
ny=tGa./hly;%isotope y
nz=tGa./hlz;%isotope z
%amount decayed (eq. 3)
dx=1./(2.^nx);%isotope x
dy=1./(2.^ny);%isotope y
dz=1./(2.^nz);%isotope z
%inital concentration (eq. 4)
C0x=x./dx;%isotope x
C0y=y./dy;%isotope z
C0z=z./dz;%isotope y
%time considered for heat production (seconds)
t=[];
%change in concentration over time considered (eq. 5)
Cxt=C0x-(C0x.*exp(-lx.*t));%isotope x
Cyt=C0y-(C0y.*exp(-ly.*t));%isotope y
Czt=C0z-(C0z.*exp(-lz.*t));%isotope z
%parameters
Na=6.02214e23;%atom/mol
Ex=[];%decay/atom for isotope x
Ey=[];%decay/atom for isotope y
Ez=[];%decay/atom for isotope z
mmx=[];%molar mass for isotope x
mmy=[];%molar mass for isotope y
mmz=[];%molar mass for isotope z
Cp=[];%specifc heat for the silicate (J/kg)
MF=[1 0.1 0.01 0.001 0.0001 0.00001];%mass fraction of layer to planet
%heat production of all isotopes considered (J) (eq. 6)
Ht=((Cxt.*Na.*Ex)./mmx)+((Cyt.*Na.*Ey)./mmy)+((Czt.*Na.*Ez)./mmz);
%temperature produced from heat production (celsius) 
Tt=(Ht./Cp)./MF;

%%
%Heat production calculation (Section 4)
%Calculation for Ht value when its equilvalent to Th (eq. 7)
Mf=[];%mass fraction of layer to planet
Hbasalt=[];%Cp*(temperature)+Lh for basalt

f=@(Ht) (Ht./Mf)-(Hbasalt.*Mf)

fzero(f,300)

%%
%Overturn modeling (Section 7)
clear all
clc
%Instability Wavelength (eq.8, lambda)
h=[];%thickness of overlying layer (m)
uunder=[];%viscosity of underlying layer (Pas)
uover=[];%viscosity of overlying layer (Pas)

lambda=2.9.*h.*(uunder./uover).^(1./3); %in m

%Instability Timescale (eq. 9, t)
rho=[];%density difference between layers (kg/m3)
gM=3.7;%Gravitational acceleration on Mercury (m/s2)

t=(6.5.*((uunder).^(2./3)).*((uover).^(1./3)))./(rho.*gM.*h); %in seconds
ty=t./(3.154.*10.^7);%conversion of time from seconds to years

%Instability Velocity (eqs. 10 and 11, a and V)
a=((3.*lambda.^2.*h)./(4.*pi)).^(1./3); %radius of instability in m
V=(1./3).*((rho.*gM.*(a.^2))./(uunder)).*((uunder+uover)./(uunder+(2./3).*uover)); %sinking velocity in m/s

%% 
%Thermal density model (Section 8)

Qcore=[0.02];%core heat flux (W/m2)
Score=[4.08.*10.^13];%surface area of core (m2)
t=[];%time in seconds
Cpdunite=[820];%specific heat of dunite (J/kgK)
r=[80];%thickness of layer being considered (m)
Mdunite=(4./3).*pi.*((2439700.^3)-((2439700-r).^3)).*3400;%mass of dunite layer (kg)
alpha=[2.*10.^-5];%thermal expansion for silicates (K-1)
rho0=[3400];%density of dunite (kg/m3)

deltarhodunite=((Qcore.*Score.*t)./(Cpdunite.*Mdunite)).*alpha.*rho0;%thermal density change of dunite (kg/m3)









