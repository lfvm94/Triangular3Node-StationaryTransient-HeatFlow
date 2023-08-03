clc
clear all

t=0.1;

%% Mesh, topology and geometry

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
alpha=10;   % W/(mK)

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

%% Heat Capacity Matrix
CapMat=zeros(ndof);
for i=1:nelements
    [ eCapMat ] = flw2ieCap( Ex(i,:), Ey(i,:), ...
                            specHeatCap(matrlIndex(i)), 4 );
    
    CapMat(Edof(i,2),Edof(i,2))=CapMat(Edof(i,2),Edof(i,2))+eCapMat(1,1);
    CapMat(Edof(i,2),Edof(i,3))=CapMat(Edof(i,2),Edof(i,3))+eCapMat(1,2);
    CapMat(Edof(i,2),Edof(i,4))=CapMat(Edof(i,2),Edof(i,4))+eCapMat(1,3);

    CapMat(Edof(i,3),Edof(i,2))=CapMat(Edof(i,3),Edof(i,2))+eCapMat(2,1);
    CapMat(Edof(i,3),Edof(i,3))=CapMat(Edof(i,3),Edof(i,3))+eCapMat(2,2);
    CapMat(Edof(i,3),Edof(i,4))=CapMat(Edof(i,3),Edof(i,4))+eCapMat(2,3);

    CapMat(Edof(i,4),Edof(i,2))=CapMat(Edof(i,4),Edof(i,2))+eCapMat(3,1);
    CapMat(Edof(i,4),Edof(i,3))=CapMat(Edof(i,4),Edof(i,3))+eCapMat(3,2);
    CapMat(Edof(i,4),Edof(i,4))=CapMat(Edof(i,4),Edof(i,4))+eCapMat(3,3);
end

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

atime=[a]; % Node Temperatures at time 0 

TempOuterBound=sum(a(boundaryEdof(1:20,2)))/length(boundaryEdof(1:20,2));

%% Transient Heat-Flow       
TempOuterTime=[TempOuterBound];
timeVector=[0];

timeSteps=2500; % total time = dt*timeSteps
fbtime=[fb];

teta=0.45;
lambda=eigen(K,CapMat);
lambda=max(lambda);
dt=abs(2/((1-2*teta)*lambda))*0.9
total_time=timeSteps*dt

% ---------------------------------------------------------------------
% New boundary conditions
% ---------------------------------------------------------------------

Tout=-(100+273.15);  %[K]

K=globalK;
fbout=zeros(ndof,1);
for i=1:20
    [ Kce,fce ] = convecte(boundaryEx(i,:), boundaryEy(i,:),alpha,t,Tout);
    [K,fbout]=assem(boundaryEdof(i,:),K,Kce,fbout,fce);

end 

fbin=zeros(ndof,1);
for i=21:143
    [ Kce, fce ] = convecte(boundaryEx(i,:), boundaryEy(i,:),alpha,t,Tin);
    [K,fbin]=assem(boundaryEdof(i,:),K,Kce,fbin,fce);

end  

A=CapMat+K*dt*teta;
fb=fbin-fbout;

% ----------------------------------------------------------------------
% Time analysis
% ----------------------------------------------------------------------

timeStepPlot=timeSteps; % Time step at which to plot the 
                        % temperatures on the wall
for i=1:timeSteps
    fbt=teta*fb+(1-teta)*fbtime(:,i);
    fbtime=[fbtime,fbt];
    
    b=(CapMat-dt*(1-teta)*K)*atime(:,i)+dt*fbt;

    an=A\b;
    atime=[atime,an];
    
    TempOuterBound=sum(an(boundaryEdof(1:20,2)))/length(boundaryEdof(1:20,2));
          
    % Warming rate of the outer temperature
    TempOuterTime=[TempOuterTime;TempOuterBound];
    timeVector=[timeVector;i];
end

%% Plot of solution
Ed=extract(Edof,atime(:,1));
figure(3)
fill(Ex',Ey',Ed')
title('Stationary 2D Wall Heat Transfer (Units: Kelvin) - Time 0 seconds')
xlabel('Width M')
ylabel('Height M')
hold on

Ed=extract(Edof,atime(:,timeStepPlot));
figure(4)
fill(Ex',Ey',Ed')
hold on
title(strcat(' Temperature of 2D Wall (Units: Kelvin) - ',...
             ' Time = ', num2str(round(dt*timeStepPlot,2)),' seconds'))
xlabel('Width M')
ylabel('Height M')
hold on

figure(5)
plot(timeVector,TempOuterTime,'b')
title('Average temperature variation rate at the outer boundary')
xlabel('Time-Steps')
ylabel('Temperature (K)')
hold on

%% Total Heat Flux
TotalHeatFluxOut=0;
elementHeatFlux=zeros(20,1);
for i=1:20
    le=((boundaryEx(i,2)-boundaryEx(i,1))^2+(boundaryEy(i,2)-...
        boundaryEy(i,1))^2)^0.5;
    elementHeatFlux(i)=alpha*t*le*((an(boundaryEdof(i,2))+...
        an(boundaryEdof(i,3)))*0.5-Tout);
    TotalHeatFluxOut=elementHeatFlux(i)+TotalHeatFluxOut;
end

disp('Total Heat-flux')
disp(TotalHeatFluxOut);

% ------------------------------ End ----------------------------------
