   % FIXED BOUNDARY PROBLEM

   % Program BD_N with Verlet neighbour lists - J19, 2012 revised
   %  Bidisperse spheres - Types A and B
   %  A has diameter eq 1 (diameter sigma = 1)  and velocity v0A
   %  VEL(N) AND RANCELL(N) - reflect this difference
   
   % Periodic boundary conditions in 2D,
   % Boxlength of {-L/2, L/2} in simulation units centered at origin 
   
   % (xpU, ypU, TpU) are unfolded coordinates
   % (xpF, ypF, TpF) are  folded coordinates
   
   % Time scaled with (3/D_r)
   
   % SET RANDOM NUMBER before using PARFOR
 
   % LOAD SAVED CONFIGURATION IF NEEDED

   % load ('./*.mat')

   % DEFINE PARAMETERS 
   
   %  cellrunID = 1;
 
   %  fileIDF = fopen('cellF1.dat','w');
   %  fileIDU = fopen('cellU1.dat','w');
alpha = 20;
createAdir = sprintf('A=%g', alpha);
mkdir(createAdir)
cd(createAdir)
for Pe = [0.1, 0.5, 1, 5, 10]
   createdir = sprintf('Pe=%g', Pe);
   mkdir(createdir)
   cd(createdir)
  [N, nu, D_R, D_T, dt, nT, sig, sizedif, r_cut, r_cutE, r_list, pre_F, pre_T, delta, Nsave] = Parameter_file(alpha);     
   
   radcell = zeros(1,N);                                                                         
   
    for Ncell = 1:N
        radcell(Ncell) = sig*0.5;
    end
       
% ALLOCATE MEMORY 
% initialize INITIAL CONFIGURATION 

    t = 0.d0;                                                              
    savecount = 0;     
%     Nhist = 1000;       % interval for count of histogram 
%     histcount = 0;      
%     nhist = nT/Nhist;
%     
%     xhist = zeros(1,nhist*N/2);     
%     yhist = zeros(1,nhist*N/2);
    
%     xhist = zeros(1,N);     
%     yhist = zeros(1,N);
    
    [xpFold, ypFold, theFold, BoxL] = diluteB(N, delta, sig);

     Boxhalf = BoxL/2.0; 
     BoxL_R = BoxL - sig;
     Boxhalf_R = BoxL_R/2;
     %      Area = BoxL*BoxL;  
     
% %  Density of cells

%     50 % is CELL A and 50% is CELL B 
%     Thus the formulae:
%     AcellA = (N/2.0)*pi*(rsigA)*(rsigA);
%     AcellB = (N/2.0)*pi*(rsigB)*(rsigB);
%     densA = AcellA/Area;
%     densB = AcellB/Area;
%     density = densA + densB;

    xp1 = xpFold;
    yp1 = ypFold;
    thep1 = theFold;
    
    tp1 = 0.0;
% INITIALIZE WALL INTERACTIONS   
    
    d1 = zeros(1,N);
    d2 = zeros(1,N);
    for cell = 1:N
       if xpFold(cell)> 0
           d1(cell) =  Boxhalf -xpFold(cell);
       else
           d1(cell) = -Boxhalf -xpFold(cell);
       end
       if ypFold(cell)> 0
           d2(cell) =  Boxhalf -ypFold(cell);
       else
           d2(cell) = -Boxhalf -ypFold(cell);
       end
    end

% CREATE LISTS FOR FIRST USE 
%     [POINT, LIST, xpFsaved, ypFsaved, theFsaved, timesaved] = LISTUB2(r_list, N, xpFold, ypFold, theFold, BoxL, t);

% INITIALIZE TOTAL FORCES AND TORQUES
  
%  % First the Linear Soring Force
          
% [LSFTx, LSFTy, LSTTz] = ForceLSB(POINT, LIST, r_cut,  N, xpFold, ypFold, BoxL, sig, dt);

           LSFTx = zeros(1,N);
           LSFTy = zeros(1,N);
           LSTTz = zeros(1,N);
 
%  % The Elastic forces
     
%   [EFTx, EFTy, ETTz] = ForceEB(POINT, LIST, r_cutE,  N, xpFold, ypFold, theFold, BoxL, pre_F, pre_T, nu, radcell);
  
           EFTx = zeros(1,N);
           EFTy = zeros(1,N);
           ETTz = zeros(1,N);

% % The Wall Interaction

%   [FxWall, FyWall, TzWall] = wall_E(d1, d2, theFold, pre_F, pre_T, N, r_cutE, nu);
  
           FxWall = zeros(1,N);
           FyWall = zeros(1,N);
           TzWall = zeros(1,N);

          % Sum to get total force

%         FTx = LSFTx + EFTx + FxWall;
%         FTy = LSFTy + EFTy + FyWall;
%         TTz = LSTTz + ETTz + TzWall;
    
 % GENERATE trial displacement for FIRST TIME STEP
  
 disX = sqrt(2.0*D_T*dt)*randn(1,N);
 disY = sqrt(2.0*D_T*dt)*randn(1,N);
 disT = sqrt(2.0*D_R*dt)*randn(1,N);
 
 % FIRST STEP - MOVE CELLS 

  xpFnew = xpFold ; 
  ypFnew = ypFold ; 
  theFnew = theFold ;
 
 
  xpU = xpFnew;         % Unfolded coordinates
  ypU = ypFnew;         % Unfolded coordinates
  theU = theFnew;       % Unfolded coordinates
 
 % FIRST STEP - ENFORCE BOUNDARY 

[xpFnew, ypFnew] = boundary(xpFnew, ypFnew, N, BoxL, sig);
     
 % UPDATE OLD VALUES FOR USE 
    
    xpFold = xpFnew;
    ypFold = ypFnew;
    theFold = theFnew;
    
% START TIME STEPPER - use ADAPTIVE EXPLICIT EULER for time stepper
tic
for nTime = 1: nT
    % UPDATE LIST if needed 
    for cell = 1:N
       if xpFold(cell)> 0
           d1(cell) =  Boxhalf -xpFold(cell);
       else
           d1(cell) = -Boxhalf -xpFold(cell);
       end
       if ypFold(cell)> 0
           d2(cell) =  Boxhalf -ypFold(cell);
       else
           d2(cell) = -Boxhalf -ypFold(cell);
       end
    end
   
  
    [POINT, LIST, xpFsaved, ypFsaved, theFsaved, timesaved] = LISTUB2(r_list, N, xpFold, ypFold, theFold, BoxL, t);
    
 % GENERATE RANDOM NUMBERS for noise 

 % DEFINE ARRAYS

  CT = cos(theFold);
  ST = sin(theFold);

  % CALCULATE TOTAL FORCE USING LISTS BASED ON OLD positions
     
   % % First the Linear Spring force

      [LSFTx, LSFTy, LSTTz] = ForceLSB(POINT, LIST, r_cut, N, xpFold, ypFold, BoxL, sig, dt);

   % % Then  Elastic forces
     
      [EFTx, EFTy, ETTz] = ForceEB(POINT, LIST, r_cutE,  N, xpFold, ypFold, theFold, BoxL, pre_F, pre_T, nu, radcell);

   % % Wall Interactions
        FrLS = 1/dt;
       [FxWall, FyWall, TzWall] = wall_E(d1, d2, FrLS, N, sig);
 
     %  SUM TO GET TOTAL FORCE and torque

        FTx = LSFTx + EFTx + FxWall;
        FTy = LSFTy + EFTy + FyWall;
        TTz = LSTTz + ETTz + TzWall;
%         FTx = 0.0;
%         FTy = 0.0;
%         TTz = 0.0;
% % DISPLACEMENT MADE 
    % Displacement due to Forces     
     FDX = dt*(FTx + Pe*CT); 
     FDY = dt*(FTy + Pe*ST); 
     TDT = dt*TTz;     
    % Displacement due to diffusion 
     disX = sqrt(2.0*D_T*dt)*randn(1,N);
     disY = sqrt(2.0*D_T*dt)*randn(1,N);   
     disT = sqrt(2.0*D_R*dt)*randn(1,N);
 
     xpU = xpU + FDX + disX; 
     ypU = ypU + FDY + disY; 
     theU = theU + TDT + disT; 
   
    xpFnew = xpFold + FDX + disX; 
    ypFnew = ypFold + FDY + disY;
    theFnew = theFold + TDT + disT;

% ENFORCE WALL/FIXED BOUNDARY     
[xpFnew, ypFnew] = boundary(xpFnew, ypFnew, N, BoxL, sig);
  % COUNT FOR HISTOGRAM
  
%   if (t>5000)
%      if (mod(nTime,Nhist) == 0)
%          for icell=1:N
%              xhist(histcount*N+icell) = xpFnew(icell);
%              yhist(histcount*N+icell) = ypFnew(icell);
%          end
%          histcount = histcount +1;
%      end
%   end
 
  % UPDATE CURRENT TIME

    t = t + dt;

  % save data, fig, mat and also binned velocity data to get averages
  % to do the binning move back to actual bins using BoxL  
                                         
 if (mod(nTime,Nsave) == 0)
   
   savecount = savecount + 1;    
   savedataB(t, xpU, ypU, theU, xpFnew, ypFnew, theFnew, Pe, D_R, D_T, LSFTx, LSFTy, LSTTz, EFTx, EFTy, ETTz, FTx, FTy, TTz, alpha, nu, N, savecount, BoxL);

  % reset values for next save and print
 
%         xp1 = xp2;
%         yp1 = yp2;
%         tp1 = tp2;

 end
         
 % UPDATE VARIABLES FOR NEXT RUN 
        
        xpFold = xpFnew;
        ypFold = ypFnew;
        theFold = theFnew;     
        
end
toc
cd ../
end
cd ../