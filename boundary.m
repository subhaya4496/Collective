function [xpFnew, ypFnew] = boundary(xpFnew, ypFnew, N, BoxL, sig)
BoxL_R = BoxL - sig;
BoxW_R = BoxL - sig;
BoxLhalf = BoxL/2;
Boxhalf_R = BoxW_R/2.0;
for icell = 1: N
% if x- direction has periodicity
  if (xpFnew(icell)  > +BoxLhalf)   
      xpFnew(icell) = xpFnew(icell) - BoxL;                
  elseif (xpFnew(icell)  < -BoxLhalf)   
      xpFnew(icell) = xpFnew(icell) + BoxL;
  end
% % Fixed Boundary in x-direction 
%   if (xpFnew(icell)  > +Boxhalf_R)   
%       xpFnew(icell) =  BoxL_R - xpFnew(icell);                
%   elseif (xpFnew(icell)  < -Boxhalf_R)   
%       xpFnew(icell) = - xpFnew(icell) - BoxL_R;                
%   end           
 
  if (ypFnew(icell)  > +Boxhalf_R)   
      ypFnew(icell) =  BoxW_R - ypFnew(icell);                
  elseif (ypFnew(icell)  < -Boxhalf_R)   
      ypFnew(icell) = - ypFnew(icell) - BoxW_R;                
  end      

end