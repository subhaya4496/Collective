% M21, 2012
 
function  [] =  savedataB(t, xpU, ypU, theU, xpFnew, ypFnew, theFnew, v0, D_R, D_T, LSFTx, LSFTy, LSTTz, EFTx, EFTy, ETTz, FTx, FTy, TTz, alpha, nu, N, savecount, BoxL)


% maindir = sprintf('V%4.2f_A%4.2f_Mu%4.2f_BR%4.2f_BT%4.2f_N%4.2f',v0,alpha, nu, B_R, B_T, N);    
  
  
%  mkdir(maindir);
% cd(maindir);
%   
      afname = sprintf('Data%1u.mat', savecount); 
   % binfname = sprintf('Bin%f.mat', savecount);
   %vorfname = sprintf('Vor%f.mat', savecount);
           
  %BDFplotname = sprintf('BDF%f.tif', savecount);
  %CDFplotname = sprintf('CDF%f.tif', savecount);
    DDFplotname = sprintf('BDM%f.tif', savecount);

  Nhalf = N/2; 
  Boxhalf = BoxL/2.0; 
  % time averaged velocities

%        XVcell = (xpFnew - xp1)./(t - tp1);
%        YVcell = (ypFnew - yp1)./(t - tp1);
 
     save(afname, 't', 'xpU', 'ypU', 'theU', 'xpFnew', 'ypFnew', 'theFnew', 'v0', 'D_R', 'D_T', 'LSFTx', 'LSFTy', 'LSTTz', 'EFTx', 'EFTy', 'ETTz', 'FTx', 'FTy', 'TTz', 'alpha', 'nu', 'N', 'savecount','BoxL');
   
   %warning('MATLAB:imagesci:tifftagsread:zeroComponentCount', 'off')
   
    figure(1)
    plot(xpFnew(1:N), ypFnew(1:N),'o', 'MarkerFaceColor','g', 'MarkerSize', 9);
    axis([ -Boxhalf  +Boxhalf  -Boxhalf  +Boxhalf]); 
    hold on;
%     plot(xpFnew(Nhalf), ypFnew(Nhalf),'o', 'MarkerFaceColor','g', 'MarkerSize', 8);
%     hold on;
%     plot(xpFnew(Nhalf+1:N-1), ypFnew(Nhalf+1:N-1),'o', 'MarkerFaceColor','y', 'MarkerSize', 8);
%     hold on;
%     plot(xpFnew(N), ypFnew(N),'o', 'MarkerFaceColor','r', 'MarkerSize', 12);
%     hold on;
    quiver(xpFnew, ypFnew, cos(theFnew), sin(theFnew), 0.1, 'LineWidth', 2, 'Color','k','ShowArrowHead','off');
    hold on
    quiver(xpFnew, ypFnew, -cos(theFnew), -sin(theFnew), 0.1, 'LineWidth', 2, 'Color','k','ShowArrowHead','off');
    hold off; 
%     quiver(xpFnew, ypFnew, XVcell, YVcell, 0); 
%     hold off;
    title(sprintf('t=%f', t));
    print(DDFplotname, '-dtiff');
   % print(DDFplotname, '-dtiff');
    %print(BDFplotname, '-dtiff');
    %saveas(gcf,DDFplotname)
    figure('visible', 'off')  
    %quiver(xpFnew, ypFnew, cos(theFnew), sin(theFnew), 0);  
    %title(sprintf('t=%f', t)); 
    %print(CDFplotname, '-dtiff');
      
 
 
    %save(vorfname, 't', 'tp1', 'v0', 'B_R', 'B_T', 'alpha', 'nu', 'N', 'savecount','BoxL', 'XVcell', 'YVcell');
    
  % quiver(xpFnew, ypFnew, XVcell, YVcell);  
  % title(sprintf('t=%f', t)); 
  % print(DDFplotname, '-dtiff');
%    saveas(gcf,DDFplotname)
  % reset values for next macroscopic velocity calculation

%     xp2 = xpFnew;
%     yp2 = ypFnew;
%     thep2 = theFnew;
% 
%     tp2 = t;
%     
    % cd;

  end

