% M30, 2012 Bidisperse

function [LSFTx, LSFTy, LSTTz] = ForceLSB(POINT, LIST, r_cut,  N, xpFold, ypFold, BoxL, sig, dt)
 
% zero forces and torques on N cells - INITIALIZE
            
           LSFTx = zeros(1,N);                                          
           LSFTy = zeros(1,N);                                         
           LSTTz = zeros(1,N);                                 
        
% Use list to find neighbours
      
for  II = 1: N-1                                                          % DO 200 I = 1, N-1
                          
                      JBEG = POINT(II);                               % corrected JBEG = POINT(I) 
                      JEND = POINT(II+1) - 1;                   % corrected JEND = POINT(I+1) -1 
          
% check that atom II has neighbours

 if (JBEG <=  JEND) 
          
                        RXI = xpFold(II);                                              
                        RYI = ypFold(II);                        
                      %  THEI = theFold(II);                                        
       
                        FXI = LSFTx(II);                                          
                        FYI = LSFTy(II);     
                        TZI = LSTTz(II);                                         
               
                      for  JNAB = JBEG: JEND                       % DO 199 JNAB = JBEG, JEND 
                    
                        JJ = LIST(JNAB);                                   % J = LIST(JNAB)
                    
                        RXIJ = RXI - xpFold(JJ);        
                        RYIJ = RYI - ypFold(JJ);                       
                      %  THEIJ = THEI - theFold(JJ);    
        
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
                  sigXcomp = sig*RXIJ/rm;
                  sigYcomp = sig*RYIJ/rm;
               % If neighbour within cut-off 
                       
               if (rm <= r_cut)                                                                          
                          
                       %     RAD_I = radcell(II);
                       %    RAD_J = radcell(JJ);
                           
                            FrLS = 1/dt;
                       % Call RLS for this pair only when rm <= r_cut   

                       FxLS = FrLS*(sigXcomp-RXIJ);
                       FyLS = FrLS*(sigYcomp-RYIJ);
                       TzLS = 0.0;                                         % LS produces no torque
                  
                       % Augment forces and torques          
                            
                       FXI = FXI + FxLS;                                   % Increment Total force x
                       FYI = FYI + FyLS;                                   % Increment Total force y
                       TZI = TZI + TzLS;                               % Increment Total torque z
                             
                       LSFTx(JJ) = LSFTx(JJ) - FxLS;                  
                       LSFTy(JJ) = LSFTy(JJ) - FyLS;                 
                       LSTTz(JJ) = LSTTz(JJ) - TzLS;                   
                  
               end
                       % end  IF cut-off loop
                       
                     end                                                 % JNAB loop (199 loop)
                                                         
                       LSFTx(II) = FXI;                                  % Total force x
                       LSFTy(II) = FYI;                                  % Total force y
                       LSTTz(II) = TZI;                                  % Total torque z 
             

 end                                                                    % END IF neighbour loop (JBEG. LE. JEND)

            
end                                                                     % end of II loop - CELL I loop

end
