function  [POINT, LIST, xpFsaved, ypFsaved, theFsaved, timesaved] = LISTUB2(r_list, N, xpFold, ypFold, theFold, BoxL, t)
                                         
           %  Point(N) Index to neighbour list 
           %  M30, 2012 Bidisperse
           %  LIST(maxnab) Verlet neighbour list
           
           % ALL DISTANCES MEASURED IN SIGMA units
           
            RLSTSQ = r_list*r_list;                                % list criterion for sorting into lists
 
           % save configuration since list is being updated
           % CALL SAVE part of this function and set this to saved output    
                            
            xpFsaved = xpFold;             
            ypFsaved = ypFold;
            theFsaved = theFold;                           
            timesaved = t;
                                 
  % START NEIGHBOUR LIST CALCULATION
           
                    NLIST = 0;                                  % reset nneighbours of cells
                         
for  CI = 1 : N-1                                            % CI = cell I   do 100 I=1, N-1           

                       POINT(CI) = NLIST + 1;                   % POINT(I) = NLIST + 1

                       RXI = xpFold(CI);                                              
                       RYI = ypFold(CI);               
                                                               
                                          
for  CJ = CI + 1 : N                                              % CJ = cell J DO 99 J=I+1, N
                                        
                       RXIJ = RXI - xpFold(CJ);        
                       RYIJ  = RYI - ypFold(CJ);
                           
   
  % MIN DISTANCE BETWEEN CI and appropriate IMAGE of CJ         

%            X distance check
 
                       if (RXIJ > +BoxL/2.0)   
                             RXIJ = RXIJ - BoxL;                
                       end 
    
                       if (RXIJ < - BoxL/2.0)   
                             RXIJ = RXIJ + BoxL;                
                       end 
              
%             % Y distance check          
%                        
%                        if (RYIJ > +BoxL/2.0)   
%                              RYIJ = RYIJ - BoxL;                
%                        end 
%     
%                        if (RYIJ < - BoxL/2.0)   
%                              RYIJ = RYIJ + BoxL;                
%                        end
                       
  % END ENFORCE MINIMUM DISTANCE
  % I just assume that MAXNAB is appropriate as in A &T terminology BASED on MaxD value
           
          RIJSQ = RXIJ*RXIJ+ RYIJ*RYIJ;                                % calculate RIJSQ in A&T terminology
                     
                      if (RIJSQ < RLSTSQ)                              % Check to see if it has to be included in LIST  
                 
                             NLIST = NLIST + 1;                         % NLIST = NLIST + 1 in A&J terminology
                             LIST(NLIST) = CJ;                          % LIST(NLIST) = J in A&J terminology
     
                      end                                               % end list cut-off loop  % endif
                                                         
end                                                                     % end cell J loop  % 99 continue 

end                                                                     %   100 continue

                  POINT(N) = NLIST + 1;                                 %  POINT(N) = NLIST + 1 in A&T terminology
                  
end                                                                     % end the function
