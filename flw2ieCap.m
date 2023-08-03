function [ elemCapMat ] = flw2ieCap( ex,ey,specHeatCap, NoIntPoints )
% function [ elemCapMat ] = flw2ieCap( ex, ey, specHeatCap,NoIntPoints ) 
%-------------------------------------------------------------------------
% Purpose: Compute the element heat capicity matrix. Supported elements:
%          a) bilinear quadrilateral and b) linear triangular
%-------------------------------------------------------------------------
% Input: ex           Element nodal x-coordinates,
%                     size(ex) = [ 1, 4 ] for bilinear quadrilaterals
%                     size(ex) = [ 1, 3 ] for linear triangles
%        ey           Element nodal y-coordinates,
%                     size(ey) = [ 1, 4 ] for bilinear quadrilaterals
%                     size(ey) = [ 1, 3 ] for linear triangles
%        specHeatCap  Specific heat capacity (rho*Cp) [J/(m^3*degC)]                     
%        NoIntPoints  Number of integration points, you may choose:
%                     1, 3 or 4 for linear triangles
%-------------------------------------------------------------------------
% Output: elemCapMat  Element heat capacity matrix
%-------------------------------------------------------------------------
% Created by: Luis F. Verduzco, 20171115
%-------------------------------------------------------------------------

% For linear triangular elements
   
% Gauss quadrature integration points and weights

if NoIntPoints == 1

    L1 = 1/3; L2 = 1/3; L3 = 1/3;
    intWeights = 1;

elseif NoIntPoints == 3

    L1 = [ 1/2 1/2 0 ]; L2 = [ 1/2 0 1/2 ]; L3 = [ 0 1/2 1/2 ];
    intWeights = [ 1/3 1/3 1/3 ];

elseif NoIntPoints == 4

    L1 = [ 1/3 0.6 0.2 0.2 ]; L2 = [ 1/3 0.2 0.6 0.2 ]; 
    L3 = [ 1/3 0.2 0.2 0.6 ];
    intWeights = [ -27/48 25/48 25/48 25/48 ];

else

     L1 = 0; L2 = 0; L3 = 0;
     intWeights = 0;
     disp( 'Used number of integration points not implemented!' );

end

intPoints = [ L1' L2' L3' ]; % Area coordinates

% Initialize element heat capacity matrix

elemCapMat = zeros( 3 );

% Loop over the integration points and compute the element heat capacity
% matrix

for indexIP = 1:NoIntPoints

    L1       = intPoints( indexIP, 1 );
    L2       = intPoints( indexIP, 2 );
    L3       = intPoints( indexIP, 3 );
    weightIP = intWeights( indexIP );

    % Write element shape functions

    Ne = [ L1 L2 L3 ];

    % Compute the Jacobian of the isopar transformation

    ex_Cross_ey = [ ex(2) - ex(1), ex(3) - ex(1);
                    ey(2) - ey(1), ey(3) - ey(1) ];                 

    trArea = 1/2*abs( det( ex_Cross_ey ) ); 

    elemCapMat = elemCapMat + specHeatCap*Ne'*Ne*trArea*weightIP;     

end

