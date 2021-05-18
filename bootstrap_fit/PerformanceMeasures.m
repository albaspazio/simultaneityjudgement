function y = PerformanceMeasures(task, LamT, LamR, Tau, Delta, Xi, PlotRange, Title)
%
% Compute manifest or latent performance measures from parameter estimates
%
% Usage: o = PerformanceMeasures(task, LamT, LamR, Tau, Delta, Xi, PlotRange, Title)
%
% Arguments:
%
% - task:      A character string. Either 'SJ2', 'SJ3', 'TOJ', or 'latent' (case
%              insensitive)
%
% - LamT:      Value of parameter lambda_t
%
% - LamR:      Value of parameter lambda_r
%
% - Tau:       Value of parameter tau
%
% - Delta:     Value of parameter delta
%
% - Xi:        Value of parameter xi
%
% - PlotRange: Option to plot the psychometric functions. A 2-element vector giving the
%              lower and upper limits of the range defining the horizontal axis for the
%              plot. An empty vector disables plotting
%
% - Title:      Plot Title. if empty, use default
%
% Output: The output is a structure with values for the applicable performance measures
%
% Further details in:  Alcalá-Quintana, R. & García-Pérez, M.A. (2013). Fitting model-based
%                         psychometric functions to simultaneity and temporal-order judgment
%                         data: MATLAB and R routines. Behavior Research Methods, in press.
%                         http://dx.doi.org/10.3758/s13428-013-0325-2

%-------------------Check conditions --------------------------------------
task = upper(task);
if sum(strcmpi(task,{'SJ2' 'SJ3' 'TOJ' 'LATENT'}))==0
  error('Invalid value in task (must be ''SJ2'' ''SJ3'' ''TOJ'' or ''latent'')')
elseif ~isscalar(LamT) || ~isnumeric(LamT) || ~isreal(LamT) || LamT<=0 
  error('Invalid value in LamT (must be a positive real scalar)')
elseif ~isscalar(LamR) || ~isnumeric(LamR) || ~isreal(LamT) || LamR<=0
  error('Invalid value in LamR (must be a positive real scalar)')
elseif ~isscalar(Tau) || ~isnumeric(Tau) || ~isreal(Tau)
  error('Invalid value in Tau (must be a real scalar)')
elseif strcmpi(task,'LATENT')==0 && ...
       (~isscalar(Delta) || ~isnumeric(Delta) || ~isreal(Delta) || Delta<0)
  error('Invalid value in Delta (must be a non-negative real scalar)')
elseif strcmpi(task,'TOJ')==1 && ...
       (~isscalar(Xi) || ~isnumeric(Xi) || ~isreal(Xi) || Xi<0 || Xi>1)
  error('Invalid value in Xi (must be a scalar within [0,1])')
elseif (numel(PlotRange)~=2 && numel(PlotRange)~=0) || ...
       (numel(PlotRange)==2 && (~isnumeric(PlotRange) || ~isreal(PlotRange) || PlotRange(1)==PlotRange(2)))
  error('Invalid PlotRange (must be empty or a two-element real vector defining a non-null range)')
elseif isempty(Title)
    Title = '\Psi_{SJ2-S}';  
end


%-------------------Main---------------------------------------------------
s1 = warning('off','MATLAB:log:logOfZero');
s2 = warning('off','MATLAB:divideByZero');
if numel(PlotRange)==2, PlotRange = sort(reshape(PlotRange,1,2)); end
tol=0.001; % accuracy for numerical computations
switch task
case 'SJ2'
  SJ2_TFsb = 'Does not exist';  SJ2_RFsb = 'Does not exist';
  SJ2_sr = 'Does not exist';  SJ2_mp = 'Does not exist';
  if Delta >= log(2)*(LamT+LamR)/(2*LamT*LamR)
    if exp(-2*Delta*LamT) <= (LamR-LamT)/(2*LamR)
      SJ2_TFsb = Delta-Tau+log((LamT+LamR)/(2*LamR*(exp(2*Delta*LamT)-1)))/LamT;
    else
      SJ2_TFsb = FindIt('TFS-2',Delta,Tau,LamT,LamR,Xi,tol);
    end
    if exp(-2*Delta*LamR) <= (LamT-LamR)/(2*LamT)
      SJ2_RFsb = log((2*LamT*(exp(2*Delta*LamR)-1))/(LamT+LamR))/LamR-Delta-Tau;
    else
      SJ2_RFsb = FindIt('RFS-2',Delta,Tau,LamT,LamR,Xi,tol);
    end
  end
  if isnumeric(SJ2_TFsb)
    SJ2_sr = SJ2_RFsb-SJ2_TFsb;
    SJ2_mp = (SJ2_RFsb+SJ2_TFsb)/2;
  end
  SJ2_peak = Delta*(LamT-LamR)/(LamT+LamR)-Tau;
  % create output structure
  y = struct('Content', 'Performance measures for SJ2 data',...
             'SJ2_TFsb', SJ2_TFsb, ...
             'SJ2_RFsb', SJ2_RFsb, ...
             'SJ2_sr', SJ2_sr, ...
             'SJ2_mp', SJ2_mp, ...
             'SJ2_peak', SJ2_peak);
case 'SJ3'
  SJ3_TFsb = 'Does not exist';  SJ3_RFsb = 'Does not exist';
  SJ3_sr = 'Does not exist';  SJ3_mp = 'Does not exist';
  if Delta >= log(2)/(2*LamR)
    if exp(-2*Delta*LamT) < (LamR-LamT)/LamT
      SJ3_TFsb = Delta-Tau+log((LamT+LamR)/(LamR*(2*exp(2*Delta*LamT)-1)))/LamT;
    else
      SJ3_TFsb = FindIt('TF-S',Delta,Tau,LamT,LamR,Xi,tol);
    end
  end
  if Delta >= log(2)/(2*LamT)
    if exp(-2*Delta*LamR) < (LamT-LamR)/LamT
      SJ3_RFsb = log((LamT*(2*exp(2*Delta*LamR)-1))/(LamT+LamR))/LamR-Delta-Tau;
    else
      SJ3_RFsb = FindIt('RF-S',Delta,Tau,LamT,LamR,Xi,tol);
    end
  end
  if isnumeric(SJ3_TFsb) && isnumeric(SJ3_RFsb)
    SJ3_sr = SJ3_RFsb-SJ3_TFsb;
    SJ3_mp = (SJ3_RFsb+SJ3_TFsb)/2;
  end
  SJ3_peak = Delta*(LamT-LamR)/(LamT+LamR)-Tau;
  if LamT <= LamR
    SJ3_TF50 = log((LamT+LamR)/(2*LamR))/LamT-Delta-Tau;
    SJ3_RF50 = Delta-Tau+log((LamT+LamR)/(2*LamR))/LamT;
  else
    SJ3_TF50 = log((2*LamT)/(LamT+LamR))/LamR-Delta-Tau;
    SJ3_RF50 = Delta-Tau-log((LamT+LamR)/(2*LamT))/LamR;
  end
  % create output structure
  y = struct('Content', 'Performance measures for SJ3 data',...
             'SJ3_TFsb', SJ3_TFsb, ...
             'SJ3_RFsb', SJ3_RFsb, ...
             'SJ3_sr', SJ3_sr, ...
             'SJ3_mp', SJ3_mp, ...
             'SJ3_peak', SJ3_peak, ...
             'SJ3_TF50', SJ3_TF50, ...
             'SJ3_RF50', SJ3_RF50);
case 'TOJ'
  ap=0.25;
  if exp(-2*LamT*Delta) >= (LamT*ap+LamR*(ap-Xi))/(LamR*(1-Xi))
    TOJ_RF25 = Delta-Tau+log(ap*(LamT+LamR)/(LamR*(1-Xi+Xi*exp(2*Delta*LamT))))/LamT;
  elseif exp(-2*LamR*Delta) >= (LamR*(1-ap)-LamT*(ap-Xi))/(LamT*Xi)
    TOJ_RF25 = log(LamT*(Xi+(1-Xi)*exp(2*Delta*LamR))/((1-ap)*(LamT+LamR)))/LamR - Delta - Tau;
  else
    TOJ_RF25 = FindIt(ap,Delta,Tau,LamT,LamR,Xi,tol);
  end
  if exp(-2*LamT*Delta) >= (LamT+LamR*(1-2*Xi))/(2*LamR*(1-Xi))
    TOJ_RF50 = Delta-Tau+log((LamT+LamR)/(2*LamR*(1-Xi+Xi*exp(2*Delta*LamT))))/LamT;
  elseif exp(-2*LamR*Delta) >= (LamR-LamT*(1-2*Xi))/(2*LamT*Xi)
    TOJ_RF50 = log(2*LamT*(Xi+(1-Xi)*exp(2*Delta*LamR))/(LamT+LamR))/LamR - Delta - Tau;
  else
    TOJ_RF50 = FindIt('PSS',Delta,Tau,LamT,LamR,Xi,tol);
  end
  ap=0.75;
  if exp(-2*LamT*Delta) >= (LamT*ap+LamR*(ap-Xi))/(LamR*(1-Xi))
    TOJ_RF75 = Delta-Tau+log(ap*(LamT+LamR)/(LamR*(1-Xi+Xi*exp(2*Delta*LamT))))/LamT;
  elseif exp(-2*LamR*Delta) >= (LamR*(1-ap)-LamT*(ap-Xi))/(LamT*Xi)
    TOJ_RF75 = log(LamT*(Xi+(1-Xi)*exp(2*Delta*LamR))/((1-ap)*(LamT+LamR)))/LamR - Delta - Tau;
  else
    TOJ_RF75 = FindIt(ap,Delta,Tau,LamT,LamR,Xi,tol);
  end
  TOJ_JND = TOJ_RF75-TOJ_RF50;
  TOJ_width = TOJ_RF75-TOJ_RF25;
  % create output structure
  y = struct('Content', 'Performance measures for TOJ data',...
             'TOJ_RF25', TOJ_RF25, ...
             'TOJ_RF50', TOJ_RF50, ...
             'TOJ_RF75', TOJ_RF75, ...
             'TOJ_JND', TOJ_JND, ...
             'TOJ_width', TOJ_width);
case 'LATENT'
  if LamR >= LamT
    latent_theta = log((LamT+LamR)/(2*LamR))/LamT - Tau;
  else
    latent_theta = log((2*LamT)/(LamT+LamR))/LamR - Tau;
  end
  if LamR/(LamT+LamR) >= 0.8413
    latent_1587 = log(0.1587*(LamT+LamR)/LamR)/LamT - Tau;
    latent_8413 = log(0.8413*(LamT+LamR)/LamR)/LamT - Tau;
  elseif LamR/(LamT+LamR) < 0.1587
    latent_1587 = log(LamT/(0.8413*(LamT+LamR)))/LamR - Tau;
    latent_8413 = log(LamT/(0.1587*(LamT+LamR)))/LamR - Tau;
  else
    latent_1587 = log(0.1587*(LamT+LamR)/LamR)/LamT - Tau;
    latent_8413 = log(LamT/(0.1587*(LamT+LamR)))/LamR - Tau;
  end
  latent_sigma = latent_8413 - latent_1587;
  % create output structure
  y = struct('Content', 'Latent measures',...
             'latent_theta', latent_theta, ...
             'latent_1587', latent_1587, ...
             'latent_8413', latent_8413, ...
             'latent_sigma', latent_sigma);
end
if numel(PlotRange)==2, PlotIt(task, PlotRange, y, LamT, LamR, Tau, Delta, Xi, Title), end
warning(s1); warning(s2);
%-----------Find numerical solution-----
function y = FindIt(string,d,Tau,lt,lr,Xi,tol)
if ischar(string)
  switch string
  case 'TFS-2'
    t=-d-Tau:tol:d*(lt-lr)/(lt+lr)-Tau;
  case 'RFS-2'
    t=d*(lt-lr)/(lt+lr)-Tau:tol:d-Tau;
  otherwise
    t=-d-Tau:tol:d-Tau;
  end
  switch string
  case 'TF-S'
    v=(2*lt*exp(-lr*(d+t+Tau))) - (lt+lr*(1-exp(-lt*(d-t-Tau))));
  case {'TFS-2','RFS-2'}
    v=(lt*exp(-lr*(d+t+Tau))) + (lr*exp(-lt*(d-t-Tau))) - (lt+lr)/2;
  case 'RF-S'
    v=(lt*exp(-lr*(d+t+Tau))) - (lt+lr*(1-2*exp(-lt*(d-t-Tau))));
  case 'PSS'
    v=lr*(1-Xi)*exp(-lt*(d-t-Tau)) - (lt*Xi*exp(-lr*(d+t+Tau))) - (lt+lr)*(1/2-Xi);
  end
elseif isnumeric(string)
  if string>0 && string<1
    t=-d-Tau:tol:d-Tau;
    v=lr*(1-Xi)*exp(-lt*(d-t-Tau)) - (lt*Xi*exp(-lr*(d+t+Tau))) - (lt+lr)*(string-Xi);
  end
end
[val pos] = min(abs(v));
y=t(pos);
% -----------------Psychometric Functions----------------------------------
function p = psi_SJ2_S(x, lam_T, lam_R, tau, delta)
  pTF = difcdf(-delta, x+tau, lam_T, lam_R);
  pRF = 1 - difcdf(delta, x+tau, lam_T, lam_R);
  pS = 1 - pTF - pRF;
  p = pS;
function p = psi_SJ3_TF(x, lam_T, lam_R, tau, delta)
  pTF = difcdf(-delta, x+tau, lam_T, lam_R);
  p = pTF;
function p = psi_SJ3_S(x, lam_T, lam_R, tau, delta)
  pTF = difcdf(-delta, x+tau, lam_T, lam_R);
  pRF = 1 - difcdf(delta, x+tau, lam_T, lam_R);
  pS = 1 - pTF - pRF;
  p = pS;
function p = psi_TOJ_TF(x, lam_T, lam_R, tau, delta, xi)
  pTF = difcdf(-delta, x+tau, lam_T, lam_R);
  pRF = 1 - difcdf(delta, x+tau, lam_T, lam_R);
  pS = 1 - pTF - pRF;
  p = pTF + (1-xi)*pS;
%------------------Distribution of the arrival-time difference-------------
function d = difcdf(x,shift,lam1,lam2)
  y = x-shift;
  left = min(lam1.*exp(lam2.*y)./(lam1+lam2),1);
  right = max(1-(lam2.*exp(-lam1.*y))./(lam1+lam2),0);
  d = (y<=0).*left + (y>0).*right;
% -----------------Plot Results--------------------------------------------
function PlotIt(task, PlotRange, y, LamT, LamR, Tau, Delta, Xi, Title)
scrsz = get(0,'ScreenSize');
figure('Position',[scrsz(3)/8 scrsz(4)/3 2*scrsz(3)/5 scrsz(4)/3])
lw = 3; % line width
fsa = 12; % font size for axis labels and title
FS = ['\fontsize{',num2str(fsa),'}'];
xinf = PlotRange(1);  xsup = PlotRange(2);
step = (xsup-xinf)/800;
xval = xinf:step:xsup; % values x-axis
top = 1.3;
HA = 'HorizontalAlignment'; VA = 'VerticalAlignment';
switch task
case 'SJ2'
  P_SJ2_S = psi_SJ2_S(xval, LamT, LamR, Tau, Delta);
  plot(xval,P_SJ2_S,'r','LineWidth',lw); hold on
  axis([xinf xsup 0 top])
  set (gca,'YTick',0:0.2:1)
  line ([xinf xsup],[1 1],'Color','k')
  h1 = 1+25*(top-1)/100; h2 = 1+66*(top-1)/100; inc = 5*(top-1)/100;
  line ([y.SJ2_peak y.SJ2_peak],[0 1],'Color','k')
  text (y.SJ2_peak,1,['peak: ',num2str(y.SJ2_peak)],HA,'Center',VA,'top','Color','k')
  if isnumeric(y.SJ2_TFsb)
    line ([y.SJ2_TFsb y.SJ2_TFsb],[0 1],'Color','k')
    text (y.SJ2_TFsb,1,['TFsb: ',num2str(y.SJ2_TFsb)],HA,'Right',VA,'bottom','Color','k')
  end
  if isnumeric(y.SJ2_RFsb)
    line ([y.SJ2_RFsb y.SJ2_RFsb],[0 1],'Color','k')
    text (y.SJ2_RFsb,1,['RFsb:',num2str(y.SJ2_RFsb)],HA,'Left',VA,'bottom','Color','k')
  end
  if isnumeric(y.SJ2_sr) 
    line ([y.SJ2_TFsb y.SJ2_RFsb],[h2 h2],'Color','k')
    line ([y.SJ2_TFsb y.SJ2_TFsb],[h2-inc h2+inc],'Color','k')
    line ([y.SJ2_RFsb y.SJ2_RFsb],[h2-inc h2+inc],'Color','k')
    text (y.SJ2_mp,h2,['range: ',num2str(y.SJ2_sr)],HA,'Center',VA,'bottom','Color','k')
    text (y.SJ2_mp,1,'\downarrow',HA,'Center',VA,'bottom','Color','k')
    text (y.SJ2_mp,h1,['midpoint: ',num2str(y.SJ2_mp)],HA,'Center',VA,'bottom','Color','k')
  end
  xlabel([FS ' Delay or SOA, \Delta\itt'])
  ylabel([FS ' Probability of judgment'])
  title([FS Title],'Color','r')
case 'SJ3'
  P_SJ3_TF = psi_SJ3_TF(xval, LamT, LamR, Tau, Delta);
  P_SJ3_S = psi_SJ3_S(xval, LamT, LamR, Tau, Delta);
  P_SJ3_RF = 1 - P_SJ3_TF - P_SJ3_S;
  plot(xval,P_SJ3_TF,'k','LineWidth',lw); hold on
  plot(xval,P_SJ3_S,'r','LineWidth',lw)
  plot(xval,P_SJ3_RF,'c','LineWidth',lw)
  axis([xinf xsup 0 top]);
  set (gca,'YTick',0:0.2:1)
  line ([xinf xsup],[1 1],'Color','k')
  h1 = 1+25*(top-1)/100; h2 = 1+66*(top-1)/100; inc = 5*(top-1)/100;
  line ([y.SJ3_peak y.SJ3_peak],[0 1],'Color','k')
  text (y.SJ3_peak,1,['peak: ',num2str(y.SJ3_peak)],HA,'Center',VA,'top','Color','k')
  if isnumeric(y.SJ3_TFsb)
    line ([y.SJ3_TFsb y.SJ3_TFsb],[0 1],'Color','k')
    text (y.SJ3_TFsb,1,['TFsb: ',num2str(y.SJ3_TFsb)],HA,'Right',VA,'bottom','Color','k')
  end
  if isnumeric(y.SJ3_RFsb)
    line ([y.SJ3_RFsb y.SJ3_RFsb],[0 1],'Color','k')
    text (y.SJ3_RFsb,1,['RFsb: ',num2str(y.SJ3_RFsb)],HA,'Left',VA,'bottom','Color','k')
  end
  if isnumeric(y.SJ3_sr)
    line ([y.SJ3_TFsb y.SJ3_RFsb],[h2 h2],'Color','k')
    line ([y.SJ3_TFsb y.SJ3_TFsb],[h2-inc h2+inc],'Color','k')
    line ([y.SJ3_RFsb y.SJ3_RFsb],[h2-inc h2+inc],'Color','k')
    text (y.SJ3_mp,h2,['range: ',num2str(y.SJ3_sr)],HA,'Center',VA,'bottom','Color','k')
    text (y.SJ3_mp,1,'\downarrow',HA,'Center',VA,'bottom','Color','k')
    text (y.SJ3_mp,h1,['midpoint: ',num2str(y.SJ3_mp)],HA,'Center',VA,'bottom','Color','k')
  end
  line ([y.SJ3_TF50 y.SJ3_TF50 xinf],[0 0.5 0.5],'Color','k')
  text (xinf,0.5,[' TF50: ',num2str(y.SJ3_TF50)],HA,'Left',VA,'bottom','Color','k')
  line ([y.SJ3_RF50 y.SJ3_RF50 xsup],[0 0.5 0.5],'Color','c')
  text (xsup,0.5,['RF50: ',num2str(y.SJ3_RF50),' '],HA,'Right',VA,'bottom','Color','c')
  xlabel([FS ' Delay or SOA, \Delta\itt'])
  ylabel([FS ' Probability of judgment'])
  title([FS '{\Psi_{SJ3-TF}               '...
              '\color{red}\Psi_{SJ3-S}               '...
              '\color{cyan}\Psi_{SJ3-RF}}'])
case 'TOJ'
  P_TOJ_TF = psi_TOJ_TF(xval, LamT, LamR, Tau, Delta, Xi);
  P_TOJ_RF = 1 - P_TOJ_TF;
  plot(xval,P_TOJ_RF,'c','LineWidth',lw); hold on
  axis([xinf xsup 0 top])
  set (gca,'YTick',0:0.2:1)
  line ([xinf xsup],[1 1],'Color','k')
  h1 = 1+33*(top-1)/100; h2 = 1+66*(top-1)/100; inc = 5*(top-1)/100;
  line ([y.TOJ_RF25 y.TOJ_RF25 xinf],[0 0.25 0.25],'Color','k')
  line ([y.TOJ_RF25 y.TOJ_RF25],[0.25 1],'LineStyle',':','Color','k')
  text (xinf,0.25,[' RF25: ',num2str(y.TOJ_RF25)],HA,'Left',VA,'bottom','Color','k')
  line ([y.TOJ_RF50 y.TOJ_RF50 xinf],[0 0.5 0.5],'Color','k')
  line ([y.TOJ_RF50 y.TOJ_RF50],[0.5 1],'LineStyle',':','Color','k')
  text (y.TOJ_RF50,0.5,['RF50 (PSS): ',num2str(y.TOJ_RF50),' '],HA,'Right',VA,'bottom','Color','k')
  line ([y.TOJ_RF75 y.TOJ_RF75 xsup],[0 0.75 0.75],'Color','k')
  line ([y.TOJ_RF75 y.TOJ_RF75],[0.75 1],'LineStyle',':','Color','k')
  text (xsup,0.75,['RF75: ',num2str(y.TOJ_RF75),' '],HA,'Right',VA,'top','Color','k')
  line ([y.TOJ_RF50 y.TOJ_RF75],[h1 h1],'Color','k')
  line ([y.TOJ_RF50 y.TOJ_RF50],[h1-inc h1+inc],'Color','k')
  line ([y.TOJ_RF75 y.TOJ_RF75],[h1-inc h1+inc],'Color','k')
  text (y.TOJ_RF75,h1,['  JND: ',num2str(y.TOJ_JND)],HA,'Left',VA,'middle','Color','k')
  line ([y.TOJ_RF25 y.TOJ_RF75],[h2 h2],'Color','k')
  line ([y.TOJ_RF25 y.TOJ_RF25],[h2-inc h2+inc],'Color','k')
  line ([y.TOJ_RF75 y.TOJ_RF75],[h2-inc h2+inc],'Color','k')
  text (y.TOJ_RF75,h2,['  width: ',num2str(y.TOJ_width)],HA,'Left',VA,'middle','Color','k')
  xlabel([FS ' Delay or SOA, \Delta\itt'])
  ylabel([FS ' Probability of judgment'])
  title([FS '\Psi_{TOJ-RF}'],'Color','c')
case 'LATENT'
  P_latent = 1 - psi_TOJ_TF(xval, LamT, LamR, Tau, 0, 0);
  plot(xval,P_latent,'b','LineWidth',lw); hold on
  axis([xinf xsup 0 top])
  set (gca,'YTick',0:0.2:1)
  line ([xinf xsup],[1 1],'Color','k')
  h1 = 1+33*(top-1)/100; h2 = 1+66*(top-1)/100; inc = 5*(top-1)/100;
  line ([y.latent_1587 y.latent_1587 xinf],[0 0.1587 0.1587],'Color','k')
  line ([y.latent_1587 y.latent_1587],[0.1587 1],'LineStyle',':','Color','k')
  text (xinf,0.1587,['  15.87%: ',num2str(y.latent_1587)],HA,'Left',VA,'bottom','Color','k')
  line ([y.latent_theta y.latent_theta xinf],[0 0.5 0.5],'Color','k')
  line ([y.latent_theta y.latent_theta],[0.5 1],'LineStyle',':','Color','k')
  text (y.latent_theta,0.5,['\theta (latent PSS): ',num2str(y.latent_theta),' '],HA,'Right',VA,'bottom','Color','k')
  line ([y.latent_8413 y.latent_8413 xsup],[0 0.8413 0.8413],'Color','k')
  line ([y.latent_8413 y.latent_8413],[0.8413 1],'LineStyle',':','Color','k')
  text (xsup,0.8413,['84.13%: ',num2str(y.latent_8413),'  '],HA,'Right',VA,'top','Color','k')
  line ([y.latent_1587 y.latent_8413],[h1 h1],'Color','k')
  line ([y.latent_1587 y.latent_1587],[h1-inc h1+inc],'Color','k')
  line ([y.latent_8413 y.latent_8413],[h1-inc h1+inc],'Color','k')
  text (y.latent_8413,h1,['  \sigma: ',num2str(y.latent_sigma)],HA,'Left',VA,'middle','Color','k')
  xlabel([FS ' Delay or SOA, \Delta\itt'])
  ylabel([FS ' Probability of judgment'])
  title([FS 'latent \Psi_{RF}'],'Color','b')
end
