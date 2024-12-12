% M30, 2012
% Incorporated the factor that converts elastic potential to LJ distances

function  [EFTx, EFTy, ETTz] = ForceEB(POINT, LIST, r_cutE,  N, xpFold, ypFold, theFold, BoxL, pre_F, pre_T, nu, radcell)
 
    % zero forces and torques on N cells
            
           EFTx = zeros(1,N);                                          
           EFTy = zeros(1,N);                                         
           ETTz = zeros(1,N);                                 
        
    % Use list to find neighbours
      
for  II = 1: N-1                                                          % DO 200 I = 1, N-1
                          
                      JBEG = POINT(II);                                    % corrected JBEG = POINT(I) 
                      JEND = POINT(II+1) - 1;                              % corrected JEND = POINT(I+1) -1 
          
    % check that atom II has neighbours

 if (JBEG <=  JEND) 
          
                        RXI = xpFold(II);                                              
                        RYI = ypFold(II);                        
                        THEI = theFold(II);                                        
       
                        FXI = EFTx(II);                                          
                        FYI = EFTy(II);     
                        TZI = ETTz(II);                                         
               
                      for  JNAB = JBEG: JEND                              % DO 199 JNAB = JBEG, JEND 
                    
                        JJ = LIST(JNAB);                                   % J = LIST(JNAB)
                      
                        RXIJ = RXI - xpFold(JJ);        
                        RYIJ = RYI - ypFold(JJ);    
                        
                        THEJ = theFold(JJ);    
        
    % Minimum imageof JJ is chosen 
                      
       % X distance

                     if (RXIJ > +BoxL/2.0)   
                         RXIJ = RXIJ - BoxL;                
                     end 
    
                     if (RXIJ < - BoxL/2.0)   
                         RXIJ = RXIJ + BoxL;                
                     end 

       % Y distance              
                     
%                      if (RYIJ > +BoxL/2.0)   
%                          RYIJ = RYIJ - BoxL;                
%                      end 
%     
%                      if (RYIJ < - BoxL/2.0)   
%                          RYIJ = RYIJ + BoxL;                
%                      end   

    % Distance metric  RIJSQ in A&T terminology         
             
                  RIJSQ = RXIJ*RXIJ+ RYIJ*RYIJ;       
                  rm = sqrt (RIJSQ);    
             
               % If neighbour within elastic cut-off  (note that force
               % value is not used in evaluating this criterion)
                       
                       if (rm <= r_cutE)                      
                                                 
                         RAD_I = radcell(II);
                         RAD_J = radcell(JJ);
     
                 % Rescale ALL LENGTHS ONCE THE CUT-OFF criterion IS SATISFIED
                           
%                           factorLJ2E = (1.0/sqrt(3.0));        % to make everything O(1) we evaluate in rescaled units
                          factorLJ2E = 1.0;
                           
                           rr = rm*factorLJ2E;                    % separation
                           x1o = RXI*factorLJ2E;                % separation - will need diff
                           x2o = (RXI-RXIJ)*factorLJ2E;       % separation - will need diff
                           y1o = RYI*factorLJ2E;                 % separation - will need diff
                           y2o = (RYI-RYIJ)*factorLJ2E;        % separation - will need diff
                           
                           theo1 = THEI;
                           theo2 = THEJ;
                           
                  % [Fx1E, Fy1E, Fx2E, Fy2E, T1E, T2E] = E_cellB(x1o,x2o,y1o,y2o,theo1,theo2, rr, con1,con2);    
                   [Fx1E, Fy1E, Fx2E, Fy2E, T1E, T2E] = E_cell_updatedB(x1o,x2o,y1o,y2o,theo1,theo2, rr, pre_F, pre_T, nu);    

                  % Formulate Switch function to cut off elastic
                  % force for small separations     
                  % effective HS Diameter is RAD_I + RAD_J 

                  r00 = 1.0*(RAD_I + RAD_J);
                  
%                   switchE = 0.5*( 1 + tanh((rm-r00)/0.01) );               % switch checked M30, 2012
if rm<r00
    switchE = 0  ;
else
    switchE = 1  ;
end
                  % Augment forces and torques          
                            
                       FXI = FXI - Fx1E*switchE;                           % Increment Total force x
                       FYI = FYI - Fy1E*switchE;                           % Increment Total force y
                       TZI = TZI + T1E*switchE;                            % Increment Total torque z
                             
                       EFTx(JJ) = EFTx(JJ) - Fx2E*switchE;                  
                       EFTy(JJ) = EFTy(JJ) - Fy2E*switchE;                 
                       ETTz(JJ) = ETTz(JJ) + T2E*switchE;                   
                  
                       end                                                 % end  IF cut-off loop
                       
                      end                                                  % JNAB loop (199 loop)
                                                         
                       EFTx(II) = FXI;                                     % Total force x
                       EFTy(II) = FYI;                                     % Total force y
                       ETTz(II) = TZI;                                     % Total torque z 
             

 end                                                                       % END IF neighbour loop (JBEG. LE. JEND)

            
end                                                                        % end of II loop - CELL I loop

end
