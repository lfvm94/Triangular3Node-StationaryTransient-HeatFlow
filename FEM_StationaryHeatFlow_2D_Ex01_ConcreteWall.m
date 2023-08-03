clc
clear all

t=0.1;

%% Mesh, geometry and topology
[Edof,Dof,Ex,Ey,boundaryEdof,boundaryEx,boundaryEy,matrlIndex]=...
importMesh('Edof.txt','Dof.txt','Ex.txt','Ey.txt','boundaryEdof.txt',...
           'boundaryEx.txt','boundaryEy.txt','matrlIndex.txt');

nelements=length(Edof);
ndof=length(Dof);

%% Materials

% thermal conductivities
kw=0.05;   % W/mK WOOL
kp=0.5;    % W/mK Plaster dense
kc=1.7;    % W/mK Concrete
kthermalFactors=[kc;kw;kp];

% heat transfer convection coefficient
alpha=10;   %W/(mK)

% ---------------------------------------------------------------------
% Density
% ---------------------------------------------------------------------
% Concrete
pc=2400*880;  %J/(m3*K)
% Wool
pw=1320*1360;       %J/(m3*K)
% Plaster
pp=700*1000;       %J/(m3*K)

specHeatCap=[pc;pw;pp];

%% Boundary conditions
Tout=-(5+273.15);  %[K]
Tin=30+273.15;   %[K]

%% Heat conduction matrix
globalK=zeros(ndof);
for i=1:nelements
    k=kthermalFactors(matrlIndex(i));
    D=[k 0;
        0 k];
    [Ke,fe]=flw2te(Ex(i,:),Ey(i,:),t,D,0);
    
    globalK(Edof(i,2),Edof(i,2))=globalK(Edof(i,2),Edof(i,2))+Ke(1,1);
    globalK(Edof(i,2),Edof(i,3))=globalK(Edof(i,2),Edof(i,3))+Ke(1,2);
    globalK(Edof(i,2),Edof(i,4))=globalK(Edof(i,2),Edof(i,4))+Ke(1,3);

    globalK(Edof(i,3),Edof(i,2))=globalK(Edof(i,3),Edof(i,2))+Ke(2,1);
    globalK(Edof(i,3),Edof(i,3))=globalK(Edof(i,3),Edof(i,3))+Ke(2,2);
    globalK(Edof(i,3),Edof(i,4))=globalK(Edof(i,3),Edof(i,4))+Ke(2,3);

    globalK(Edof(i,4),Edof(i,2))=globalK(Edof(i,4),Edof(i,2))+Ke(3,1);
    globalK(Edof(i,4),Edof(i,3))=globalK(Edof(i,4),Edof(i,3))+Ke(3,2);
    globalK(Edof(i,4),Edof(i,4))=globalK(Edof(i,4),Edof(i,4))+Ke(3,3);

end

%%                      Heat convection matrix 
%%%%%%%%%%%%%%-----------------Kc & fbc----------------%%%%%%%%%%%%%%%%%%
%%%-------------------------Inner Boundary---------------------------%%%
K=zeros(ndof);
fbin=zeros(ndof,1);
for i=21:143
    [ Kce, fce ] = convecte(boundaryEx(i,:),boundaryEy(i,:),alpha,t,Tin );
    [K,fbin]=assem(boundaryEdof(i,:),K,Kce,fbin,fce);

end  
%-----------------------------------------------------------------------
K=globalK+K;
%%%-------------------------Outter Boundary---------------------------%%%
fbout=zeros(ndof,1);
for i=1:20
    [ Kce,fce ] = convecte(boundaryEx(i,:),boundaryEy(i,:),alpha,t,Tout);
    [K,fbout]=assem(boundaryEdof(i,:),K,Kce,fbout,fce);

end 
fb=fbin-fbout;

%% Solving the system of equations
[a,r]=solveq(K,fb);

%% Plot of solution
Ed=extract(Edof,a);
figure(3)
fill(Ex',Ey',Ed')
title('Stationary 2D Wall Heat Transfer (Units: Kelvin)')
xlabel('Width M')
ylabel('Height M')
hold on

%% Total Heat Flux
TotalHeatFluxOut=0;
elementHeatFlux=zeros(20,1);
TempsLeftBoundary=zeros(20,2);
for i=1:20
    le=((boundaryEx(i,2)-boundaryEx(i,1))^2+(boundaryEy(i,2)-boundaryEy(i,1))^2)^0.5;
    elementHeatFlux(i)=alpha*t*le*((a(boundaryEdof(i,2))+a(boundaryEdof(i,3)))*0.5-Tout);
    TempsLeftBoundary(i,1)=a(boundaryEdof(i,2));
    TempsLeftBoundary(i,2)=a(boundaryEdof(i,3));
    TotalHeatFluxOut=elementHeatFlux(i)+TotalHeatFluxOut;
end

disp('Total Heat-flux')
disp(TotalHeatFluxOut);

% -------------------------------- End ----------------------------------