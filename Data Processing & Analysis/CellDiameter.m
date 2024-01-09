function [cellDiameters] = CellDiameter(deltaR, R, Deff, L)
% Determines cell diameter [d] using the equation:
%
%  deltaR       d^3                   1
%  ------ = ------------ * [ ------------------ ]
%    R       (Deff^2*L)      1 - 0.8*(d/Deff)^3
%
% Equation is from O. A. Saleh's Thesis (2003), specifically for spherical
% particles of intermediate size range relative to the effective pore
% diameter. This equation is subsequently utlized by J. H. Kim in
% "Mechano-node-pore Sensing" (2018) and "Visco-node-pore Sensing" (2019).
%  
% Input:
% [deltaR] : Vector of measured delta R in pore or contraction [Ohms]
% [R] : Vector of baseline resistance (in Zone #) [Ohms]
% [Deff] : Vector of effective diameters of pore or contraction channel [um]
% [L] : Total length of channel [um]
%
% Output:
% [cellDiameters] : a vector of calculated cell diameters in [um]

% initialize output vector
cellDiameters = zeros(length(deltaR),1); 

% 
for i = 1:length(deltaR)
    
    cellDiameters(i) = ( (deltaR(i) * (Deff(i)^2) * L ) / ...
                         ( ((deltaR(i) * 0.8 * L) / Deff(i)) + R(i) ) ...
                       ) ^(1/3);
                                          
end

