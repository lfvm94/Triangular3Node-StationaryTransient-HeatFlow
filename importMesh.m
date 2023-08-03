function [Edof,Dof,Ex,Ey,boundaryEdof,boundaryEx,boundaryEy,matrlIndex]=...
         importMesh(Edof,Dof,Ex,Ey,boundaryEdof,boundaryEx,boundaryEy,...
         matrlIndex)
%-------------------------------------------------------------------------
% function [Edof,Dof,Ex,Ey,boundaryEdof,boundaryEx,boundaryEy,matrlIndex]=...
%          importMesh(Edof,Dof,Ex,Ey,boundaryEdof,boundaryEx,boundaryEy,...
%          matrlIndex) 
%-------------------------------------------------------------------------
% Purpose: To import the mesh data from .txt files
%-------------------------------------------------------------------------
% Input: Edof           name of file containing the topology matrx
%                       'Edof.txt'
%
%
%
%-------------------------------------------------------------------------
% Output: elemCapMat  Element heat capacity matrix
%-------------------------------------------------------------------------
% Created by: Luis F. Verduzco, 20171115
%-------------------------------------------------------------------------
%% Importing mesh data
nombre_archivo=Edof;
fid=fopen(nombre_archivo);
Edof=importdata(nombre_archivo);

nombre_archivo=Dof;
fid=fopen(nombre_archivo);
Dof=importdata(nombre_archivo);

nombre_archivo=Ex;
fid=fopen(nombre_archivo);
Ex=importdata(nombre_archivo);

nombre_archivo=Ey;
fid=fopen(nombre_archivo);
Ey=importdata(nombre_archivo);

nombre_archivo=boundaryEdof;
fid=fopen(nombre_archivo);
boundaryEdof=importdata(nombre_archivo);

nombre_archivo=boundaryEx;
fid=fopen(nombre_archivo);
boundaryEx=importdata(nombre_archivo);

nombre_archivo=boundaryEy;
fid=fopen(nombre_archivo);
boundaryEy=importdata(nombre_archivo);

nombre_archivo=matrlIndex;
fid=fopen(nombre_archivo);
matrlIndex=importdata(nombre_archivo);

% ----------------------------- End --------------------------------------