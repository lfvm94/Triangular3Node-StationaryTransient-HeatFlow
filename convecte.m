function [ Kce, fce ] = convecte( ex, ey, alpha, thick, Tamb )
% function [ Kce, fce ] = convecte( ex, ey, alpha, thick, Tamb )
%-------------------------------------------------------------------
% Purpose: Compute the element stiffness matrix and force vector
%          due to convection for a linear boundary segment
%-------------------------------------------------------------------
% Input: ex      Element nodal x-coordinate, size(ex) = [1, 2]
%        ey      Element nodal y-coordinate, size(ey) = [1, 2]
%        alpha   Coefficient of thermal expansion, [W/(m^2*degC)]
%        thick   Out-of-plane thickness of the continuum
%        Tamb    Ambient temperature at the boundary segment, [degC]
%-------------------------------------------------------------------
% Output: Kce    Element stiffness matrix due to convection
%         fce    Element force vector due to convection
%-------------------------------------------------------------------
% Created by: Luis F. Verduzco
%-------------------------------------------------------------------

lengthBoundary=((ex(2)-ex(1))^2+...
        (ey(2)-ey(1))^2)^0.5;

Kce=[alpha*thick*lengthBoundary,0;
    0,alpha*thick*lengthBoundary];

fce=[alpha*thick*Tamb*lengthBoundary;
    alpha*thick*Tamb*lengthBoundary];



