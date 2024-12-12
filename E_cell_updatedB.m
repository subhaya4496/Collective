function  [Fx1E, Fy1E, Fx2E, Fy2E, T1E, T2E] = E_cell_updatedB(x1o, x2o, y1o,y2o,theo1,theo2, rr, pre_F,pre_T,nu) 

    dx21 = x2o - x1o;
    dy21 = y2o - y1o;
    r = rr;
    R4 = (1./r)^4;
    R3 = (1./r)^3;
    alp = theo1 - theo2;
%     if (dx21 >= 0)
%        beta = asin(dy21/r);
%     elseif (dy21>=0)
%        beta = asin(dy21/r)+pi;
%     else
%        beta = pi-asin(dy21/r);
%     end
    beta = atan2(dy21, dx21);
    the = theo1 - beta;
    thep = theo2 - beta;
    sumthe = the+ thep;
    F_E = (pre_F*R4)*(2*(nu-1.)+6*(nu-1.)*cos(2.*the)+(nu-2.)*cos(2.*alp)+6.*(nu-1.)*cos(2.*thep)-15.*nu*cos(2.*sumthe));
  % F2E = -F1E;
    Fx1E = F_E*cos(beta) - (pre_T*R4)*(6*(nu-1)*(sin(2*the) + sin(2*thep))-30*nu*sin(2*sumthe))*sin(beta);
    Fx2E = -Fx1E;
    Fy1E = F_E*sin(beta) + (pre_T*R4)*(6*(nu-1)*(sin(2*the) + sin(2*thep))-30*nu*sin(2*sumthe))*cos(beta);
    Fy2E = -Fy1E;
    T1E = (pre_T*R3)*(-6.*(nu-1.)*sin(2.*the)-(nu-2.)*sin(2.*alp)+15.*nu*sin(2.*sumthe));
    T2E = (pre_T*R3)*(-6.*(nu-1.)*sin(2.*thep)+(nu-2.)*sin(2.*alp)+15.*nu*sin(2.*sumthe));
end