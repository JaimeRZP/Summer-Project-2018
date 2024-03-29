function [f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2,ss]= ...
      leasqr(x,y,pin,F,stol,niter,wt,dp,dFdp,options)
%function [f,p,kvg,iter,corp,covp,covr,stdresid,Z,r2,ss]=
%                   leasqr(x,y,pin,F,{stol,niter,wt,dp,dFdp,options})
%
% Levenberg-Marquardt nonlinear regression of f(x,p) to y(x).
%
% Version 3.beta
% Optional parameters are in braces {}.
% x = column vector or matrix of independent variables, 1 row per
%   observation: x = [x0 x1....xm].
% y = column vector of observed values, same number of rows as x.
% wt = column vector (dim=length(x)) of statistical weights.  These
%   should be set to be proportional to (sqrt of var(y))^-1; (That is,
%   the covariance matrix of the data is assumed to be proportional to
%   diagonal with diagonal equal to (wt.^2)^-1.  The constant of
%   proportionality will be estimated.); default = ones(length(y),1).
% pin = column vec of initial parameters to be adjusted by leasqr.
% dp = fractional increment of p for numerical partial derivatives;
%   default = .001*ones(size(pin))
%   dp(j) > 0 means central differences on j-th parameter p(j).
%   dp(j) < 0 means one-sided differences on j-th parameter p(j).
%   dp(j) = 0 holds p(j) fixed i.e. leasqr wont change initial guess: pin(j)
% F = name of function in quotes; the function shall be of the form y=f(x,p),
%   with y, x, p of the form y, x, pin as described above.
% dFdp = name of partial derivative function in quotes; default is "dfdp", a
%   slow but general partial derivatives function; the function shall be
%   of the form prt=dfdp(x,f,p,dp,F) (see dfdp.m).
% stol = scalar tolerance on fractional improvement in scalar sum of
%   squares = sum((wt.*(y-f))^2); default stol = .0001;
% niter = scalar maximum number of iterations; default = 20;
% options = matrix of n rows (same number of rows as pin) containing
%   column 1: desired fractional precision in parameter estimates.
%     Iterations are terminated if change in parameter vector (chg) on two
%     consecutive iterations is less than their corresponding elements
%     in options(:,1).  [ie. all(abs(chg*current parm est) < options(:,1))
%      on two consecutive iterations.], default = zeros().
%   column 2: maximum fractional step change in parameter vector.
%     Fractional change in elements of parameter vector is constrained to be 
%     at most options(:,2) between sucessive iterations.
%     [ie. abs(chg(i))=abs(min([chg(i) options(i,2)*current param estimate])).],
%     default = Inf*ones().
%
%          OUTPUT VARIABLES
% f = column vector of values computed: f = F(x,p).
% p = column vector trial or final parameters. i.e, the solution.
% kvg = scalar: = 1 if convergence, = 0 otherwise.
% iter = scalar number of iterations used.
% corp = correlation matrix for parameters.
% covp = covariance matrix of the parameters.
% covr = diag(covariance matrix of the residuals).
% stdresid = standardized residuals.
% Z = matrix that defines confidence region (see comments in the source).
% r2 = coefficient of multiple determination.
% ss=scalar sum of squares=sum-over-i(wt(i)*(y(i)-f(i)))^2.
%
% All Zero guesses not acceptable

% A modified version of Levenberg-Marquardt
% Non-Linear Regression program previously submitted by R.Schrager.
% This version corrects an error in that version and also provides
% an easier to use version with automatic numerical calculation of
% the Jacobian Matrix. In addition, this version calculates statistics
% such as correlation, etc....
%
% Version 3 Notes
% Errors in the original version submitted by Shrager (now called version 1)
% and the improved version of Jutan (now called version 2) have been corrected.
% Additional features, statistical tests, and documentation have also been
% included along with an example of usage.  BEWARE: Some the the input and
% output arguments were changed from the previous version.
%
%     Ray Muzic     <rfm2@ds2.uh.cwru.edu>
%     Arthur Jutan  <jutan@charon.engga.uwo.ca>

% Richard I. Shrager (301)-496-1122
% Modified by A.Jutan (519)-679-2111
% Modified by Ray Muzic 14-Jul-1992
%       1) add maxstep feature for limiting changes in parameter estimates
%          at each step.
%       2) remove forced columnization of x (x=x(:)) at beginning. x could be
%          a matrix with the ith row of containing values of the 
%          independent variables at the ith observation.
%       3) add verbose option
%       4) add optional return arguments covp, stdresid, chi2
%       5) revise estimates of corp, stdev
% Modified by Ray Muzic 11-Oct-1992
%	1) revise estimate of Vy.  remove chi2, add Z as return values
% Modified by Ray Muzic 7-Jan-1994
%       1) Replace ones(x) with a construct that is compatible with versions
%          newer and older than v 4.1.
%       2) Added global declaration of verbose (needed for newer than v4.x)
%       3) Replace return value var, the variance of the residuals with covr,
%          the covariance matrix of the residuals.
%       4) Introduce options as 10th input argument.  Include
%          convergence criteria and maxstep in it.
%       5) Correct calculation of xtx which affects coveraince estimate.
%       6) Eliminate stdev (estimate of standard deviation of parameter
%          estimates) from the return values.  The covp is a much more
%          meaningful expression of precision because it specifies a confidence
%          region in contrast to a confidence interval..  If needed, however,
%          stdev may be calculated as stdev=sqrt(diag(covp)).
%       7) Change the order of the return values to a more logical order.
%       8) Change to more efficent algorithm of Bard for selecting epsL.
%       9) Tighten up memory usage by making use of sparse matrices (if 
%          MATLAB version >= 4.0) in computation of covp, corp, stdresid.
% Modified by Sean Brennan 17-May-1994
%          verbose is now a vector:
%          verbose(1) controls output of results
%          verbose(2) controls plotting intermediate results
% Modified by Francesco Potortì
%       for use in Octave
%
% References:
% Bard, Nonlinear Parameter Estimation, Academic Press, 1974.
% Draper and Smith, Applied Regression Analysis, John Wiley and Sons, 1981.
%
%set default args

% argument processing
%

%if (sscanf(version,'%f') >= 4),
vernum= sscanf(version,'%f');
if vernum(1) >= 4,
  global verbose
  plotcmd='plot(x(:,1),y,''+'',x(:,1),f); figure(gcf)';
else
  plotcmd='plot(x(:,1),y,''+'',x(:,1),f); shg';
end;
if (exist('OCTAVE_VERSION'))
  global verbose
  plotcmd='plot(x(:,1),y,"+;data;",x(:,1),f,";fit;");';
end;

if(exist('verbose')~=1), %If verbose undefined, print nothing
	verbose(1)=0;    %This will not tell them the results
	verbose(2)=0;    %This will not replot each loop
end;

if (nargin <= 8), dFdp='dfdp'; end;
if (nargin <= 7), dp=.001*(pin*0+1); end; %DT
if (nargin <= 6), wt=ones(length(y),1); end;	% SMB modification
if (nargin <= 5), niter=20; end;
if (nargin == 4), stol=.0001; end;
%

y=y(:); wt=wt(:); pin=pin(:); dp=dp(:); %change all vectors to columns
% check data vectors- same length?
m=length(y); n=length(pin); p=pin;[m1,m2]=size(x);
if m1~=m ,error('input(x)/output(y) data must have same number of rows ') ,end;

if (nargin <= 9), 
  options=[zeros(n,1), Inf*ones(n,1)];
  nor = n; noc = 2;
else
  [nor, noc]=size(options);
  if (nor ~= n),
    error('options and parameter matrices must have same number of rows'),
  end;
  if (noc ~= 2),
    options=[options(:,1), Inf*ones(nor,1)];
  end;
end;
pprec=options(:,1);
maxstep=options(:,2);
%

% set up for iterations
%
f=feval(F,x,p); fbest=f; pbest=p;
r=wt.*(y-f);
ss=r'*r;
sbest=ss;
nrm=zeros(n,1);
chgprev=Inf*ones(n,1);
kvg=0;
epsLlast=1;
epstab=[.1, 1, 1e2, 1e4, 1e6];

% do iterations
%
for iter=1:niter,
  pprev=pbest;
  prt=feval(dFdp,x,fbest,pprev,dp,F);
  r=wt.*(y-fbest);
  sprev=sbest;
  sgoal=(1-stol)*sprev;
  for j=1:n,
    if dp(j)==0,
      nrm(j)=0;
    else
      prt(:,j)=wt.*prt(:,j);
      nrm(j)=prt(:,j)'*prt(:,j);
      if nrm(j)>0,
        nrm(j)=1/sqrt(nrm(j));
      end;
    end
    prt(:,j)=nrm(j)*prt(:,j);
  end;
% above loop could ? be replaced by:
% prt=prt.*wt(:,ones(1,n)); 
% nrm=dp./sqrt(diag(prt'*prt)); 
% prt=prt.*nrm(:,ones(1,m))';
  [prt,s,v]=svd(prt,0);
  s=diag(s);
  g=prt'*r;
  for jjj=1:length(epstab),
    epsL = max(epsLlast*epstab(jjj),1e-7);
    se=sqrt((s.*s)+epsL);
    gse=g./se;
    chg=((v*gse).*nrm);
%   check the change constraints and apply as necessary
    ochg=chg;
    idx = ~isinf(maxstep);
    limit = abs(maxstep(idx).*pprev(idx));
    chg(idx) = min(max(chg(idx),-limit),limit);
    if (verbose(1) & any(ochg ~= chg)),
      disp(['Change in parameter(s): ', ...
         sprintf('%d ',find(ochg ~= chg)), 'were constrained']);
    end;
    aprec=abs(pprec.*pbest);       %---
% ss=scalar sum of squares=sum((wt.*(y-f))^2).
    if (any(abs(chg) > 0.1*aprec)),%---  % only worth evaluating function if
      p=chg+pprev;                       % there is some non-miniscule change
      f=feval(F,x,p);
      r=wt.*(y-f);
      ss=r'*r;
      if ss<sbest,
        pbest=p;
        fbest=f;
        sbest=ss;
      end;
      if ss<=sgoal,
        break;
      end;
    end;                          %---
  end;
  epsLlast = epsL;
  if (verbose(2)),
    eval(plotcmd);
  end;
  if ss<eps,
    break;
  end
  aprec=abs(pprec.*pbest);
%  [aprec, chg, chgprev]
  if (all(abs(chg) < aprec) & all(abs(chgprev) < aprec)),
    kvg=1;
    if (verbose(1)),
      fprintf('Parameter changes converged to specified precision\n');
    end;
    break;
  else
    chgprev=chg;
  end;
  if ss>sgoal,
    break;
  end;
end;

% set return values
%
p=pbest;
f=fbest;
ss=sbest;
kvg=((sbest>sgoal)|(sbest<=eps)|kvg);
if kvg ~= 1 , disp(' CONVERGENCE NOT ACHIEVED! '), end;

% CALC VARIANCE COV MATRIX AND CORRELATION MATRIX OF PARAMETERS
% re-evaluate the Jacobian at optimal values
jac=feval(dFdp,x,f,p,dp,F);
msk = dp ~= 0;
n = sum(msk);           % reduce n to equal number of estimated parameters
jac = jac(:, msk);	% use only fitted parameters

%% following section is Ray Muzic's estimate for covariance and correlation
%% assuming covariance of data is a diagonal matrix proportional to
%% diag(1/wt.^2).  
%% cov matrix of data est. from Bard Eq. 7-5-13, and Row 1 Table 5.1 

if exist('sparse')  % save memory
  Q=sparse(1:m,1:m,1./wt.^2);
  Qinv=sparse(1:m,1:m,wt.^2);
else
  Q=diag((0*wt+1)./(wt.^2));
  Qinv=diag(wt.*wt);
end
resid=y-f;                                    %un-weighted residuals
covr=resid'*Qinv*resid*Q/(m-n);                 %covariance of residuals
Vy=1/(1-n/m)*covr;  % Eq. 7-13-22, Bard         %covariance of the data 

jtgjinv=inv(jac'*Qinv*jac);			%argument of inv may be singular
covp=jtgjinv*jac'*Qinv*Vy*Qinv*jac*jtgjinv; % Eq. 7-5-13, Bard %cov of parm est
d=sqrt(abs(diag(covp)));
corp=covp./(d*d');

if exist('sparse')
  covr=spdiags(covr,0);
  stdresid=resid./sqrt(spdiags(Vy,0));
else
  covr=diag(covr);                 % convert returned values to compact storage
  stdresid=resid./sqrt(diag(Vy));  % compute then convert for compact storage
end
Z=((m-n)*jac'*Qinv*jac)/(n*resid'*Qinv*resid);

%%% alt. est. of cov. mat. of parm.:(Delforge, Circulation, 82:1494-1504, 1990
%%disp('Alternate estimate of cov. of param. est.')
%%acovp=resid'*Qinv*resid/(m-n)*jtgjinv

%Calculate R^2 (Ref Draper & Smith p.46)
%
r=corrcoef([y(:),f(:)]);
r2=r(1,2).^2;

% if someone has asked for it, let them have it
%
if (verbose(2)), eval(plotcmd); end,
if (verbose(1)),
  disp(' Least Squares Estimates of Parameters')
  disp(p')
  disp(' Correlation matrix of parameters estimated')
  disp(corp)
  disp(' Covariance matrix of parameters estimated')
  disp(covp)
%  disp(' Covariance matrix of Residuals' )
%  disp(covr)
  disp(' Correlation Coefficient R^2')
  disp(r2)
  disp(['95% conf region: F(0.05)' ...
       sprintf('(%.0f,%.0f)>= delta_pvec''*Z*delta_pvec',n,m-n)]);
  disp(Z)
%   runs test according to Bard. p 201.
  n1 = sum((f-y) < 0);
  n2 = sum((f-y) > 0);
  nrun=sum(abs(diff((f-y)<0)))+1;
  if ((n1>10)&(n2>10)), % sufficent data for test?
    zed=(nrun-(2*n1*n2/(n1+n2)+1)+0.5)/(2*n1*n2*(2*n1*n2-n1-n2)...
      /((n1+n2)^2*(n1+n2-1)));
    if (zed < 0),
      prob = erfc(-zed/sqrt(2))/2*100;
      disp([num2str(prob),'% chance of fewer than ',num2str(nrun),' runs.']);
    else,
      prob = erfc(zed/sqrt(2))/2*100;
      disp([num2str(prob),'% chance of greater than ',num2str(nrun),' runs.']);
    end;
  end;
end;
