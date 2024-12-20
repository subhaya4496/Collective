function [xpFold, ypFold, theFold, BoxL] = diluteB(N, delta, sig)

          % Initial set-up for DILUTE BIDERSPERSE systems - M30, 2012                                                               
          % Calculate box size and density in SIGMA units 
          % Separation is center to center
          % delta = surface to surface separation (INITIAL)
          % all lengths are in SIGMA units
%     dis = delta/2.0; % This is arbitrary cell to cell surface separation  

%   Initialize memory and values

    xsigma = zeros(1,N);
    ysigma = zeros(1,N);
    
    xpFold = zeros(1,N);
    ypFold = zeros(1,N);
    theFold= zeros(1,N);
   %   Box geometry

   BoxL = 31.0*(sig+delta);     % BoxL chosen in SIGMA units                     
%    BoxL = 10.0*(sig+delta);
    Boxhalf = BoxL/2.0;                                     

% %  Density - DILUTE 
%   corrected for the bidisperse system

% [-BoxL/2, +BoxL/2.0] centered at the origin

 %Assign random angular orientations [0, 2 pi]
%       theFold = [0.0 pi/2.0];
    for  ncell = 1: N
          theFold(ncell) = 2*pi*rand;
    end
%Random position of cells  
    for cellcount = 1: N
        xsigma(cellcount) = (BoxL - sig)*rand;               
        ysigma(cellcount) = (BoxL - sig)*rand;
        for i = 1: cellcount
            if (abs(xsigma(cellcount)- xsigma(i))< sig) && (abs(ysigma(cellcount)- ysigma(i))< sig)
                break;
            end
        end
        cellcount = cellcount - 1;
    end

% %Start placing cells on regular lattice / grid - DILUTE SYSTEMS ONLY
%   cell index is based on (P,Q)
%   cellcount = 1;
%    for P = 1: N  
%           for Q = 1: M;
%                 xsigma(cellcount) = (P+0.5)* BoxL/M;               
%                 ysigma(cellcount) = (Q+0.5)* BoxL/M;         
%                 cellcount = cellcount + 1;
%           end
%    end
   
% Generate folded box-centric co-ordinates for simulation
%    xpFold = [0.0 1.0];
%    ypFold = [0.0 2.0];
  for cell = 1: N
   xpFold(cell) = real((xsigma(cell) - Boxhalf));
   ypFold(cell) = real((ysigma(cell) - Boxhalf));
  end
 
end