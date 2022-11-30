# Robust Estimation and Inference for High-Dimensional IV Cointegration Models 

### Description 

The R package [‘Robust_high_dimensional_IV_Cointegration’](https://github.com/christiskatsouris/Robust_high_dimensional_IV_Cointegration) (under development by Christis G. Katsouris) implements robust econometric estimation and inference methodologies for high-Dimensional IV Cointegration Models with autoregressive roots under both the regimes of stationarity and nonstationarity and persistence types as defined by Phillips and Magdalinos (2009): "Econometric inference in the vicinity of unity" (see, also Magdalinos and Phillips (2020)). The current package builds on the [‘ivxPredictive’](https://github.com/christiskatsouris/ivxPredictive) package prepared by Christis G. Katsouris. 

<p align="center">
  
<img src="https://github.com/christiskatsouris/ivxPredictive/blob/main/data/persistence.jpg" width="460"/>

</p>  

We consider the following persistence classes:
   
#### (P1). nearly stable processes: 

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \zeta = - \infty \ \ \text{and it holds that} \ \ \theta_n \to | \theta | < 1.$$
    
#### (P2). nearly unstable processes:   

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \ \zeta \equiv c \in \mathbb{R} \ \ \text{and it holds that} \ \ \theta_n \to \theta = 1.$$

    
#### (P3). nearly explosive processes:   

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \ \zeta = + \infty \ \ \text{and it holds that} \ \ \theta_n \to | \theta | > 1.$$

  
### Methodology  
  
The [‘ivxPredictive’](https://github.com/christiskatsouris/ivxPredictive) R package implements a novel endogenous instrumentation approach based on the IVX estimator proposed by [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) and [Kostakis, Magdalinos and Stamatogiannis (2015)](https://academic.oup.com/rfs/article/28/5/1506/1867633?login=true). Notice that the current procedure has a similar construction to the IV instrumentation proposed in the recent working paper of Magdalinos and Petrova (2022), with the aim to provide uniform and robust to the nuisance parameter of persistence inference across the spectrum of stationary and nonstationary roots, specifically for quantile autoregressive processes. We call this variant of the original IVX estimator, IVX-P, which can be employed to both conditional mean and conditional quantile functional forms when the model includes either univariate or multivariate regressors. The novelty of the IVX-P estimator is that is a 'hybrid estimator' which in contrast to the classical least squares estimator has desirable asymptotic theory properties and is constructed based on the underline nonstationary stochastic processes using information both within the admissible parameter space as well as outside the usual parameter space of autoregressive prcosesses.  

Furthermore, this R package implements robust estimation and testing for high-Dimensional IV Cointegration Models with either a conditional mean or conditional quantile specification function. Specifically, the user can choose the type of the functional form of the cointegration model. 
  
## Installation (under development) 

The R package [‘Robust_high_dimensional_IV_Cointegration’](https://github.com/christiskatsouris/Robust_high_dimensional_IV_Cointegration) will be able to be installed from Github.

## Usage: 

```R

# After development the package will be able to be installed using
install.packages("Robust_high_dimensional_IV_Cointegration")
library("Robust_high_dimensional_IV_Cointegration")

```

## Estimation and Inference Examples:

An important application of the proposed estimation and inference framework is the monitoring of financial bubbles. However, for the empirical application of the paper we focus on a climate change related dataset. We begin by considering the estimation and testing procedure for the proposed modelling framework. 

```MATLAB

% Procedure borrowed from the paper: Ke-Li Xu & Junjie Guo (2022). 
% "A New Test for Multiple Predictive Regression". Journal of Financial Econometrics.

function [aa bb cc dd ee ff ] = compute_size_0608(dataset,nr,h,epsilon);

%start simulation
Ys      = dataset.Y;
Xs      = dataset.X;
nr      = size(Ys,3);
pv      = zeros(nr, 1); %IVX-MAX0.

rej     = zeros(nr, 1);
rej_OLS = zeros(nr, 1);
rej_KMS = zeros(nr, 1);
rej_LM  = zeros(nr, 1);
rej_W   = zeros(nr, 1);
pv_W0   = zeros(nr, 1);
pv_W1   = zeros(nr, 1);
parfor i = 1:nr
    yt = squeeze(Ys(:,1,i));
    xt = squeeze(Xs(:,:,i));
    [rej(i,1) pv(i,1) rej_OLS(i,1) rej_W(i,1) pv_W0(i,1) pv_W1(i,1)] =  IV_max_J_IVX_0404(yt,xt,h,epsilon);
    [~,~,Wivx,~,~,~]=ivxlh_0404(yt,xt,h);
    rej_KMS(i,1) = Wivx(2);
    [~,~,Wivx,~,~,~]=ivxlh_modified_0419(yt,xt,h,0);
    rej_LM(i,1) = Wivx(2);    
end


aa=(sum(rej_KMS<=epsilon)/nr)';
bb=(sum(rej_LM<=epsilon)/nr)';
cc=(sum(pv<=epsilon)/nr);
dd=(sum(rej)/nr);
ee=(sum(pv_W0<=epsilon)/nr);
ff=(sum(pv_W1<=epsilon)/nr);


end

function [rej,pv,rej1,rej_W,pv_W0, pv_W1] = IV_max_J_IVX_0404(R,X,h,epsilon)
%% housekeeping
T               = size(X,1);
m               = size(X,2);
mf              = 0;
B               = 10000;

%% start OLS regression
YY              = R(2:end);
XX              = [X(1:end-1,:)];
nn              = size(YY,1);
YY              = YY - mean(YY);
XX              = XX - kron(ones(nn,1),mean(XX));
Yvec            = kron(YY,ones(m,1));

%construct Xmat
Xmat            = zeros(nn*m,m*(1+mf));
for i = 1:nn
    for j = 1:m
        Xmat((i-1)*m+j,1+(j-1)*(1+mf))  = [XX(i,j)];
    end
end
theta           = (Xmat'*Xmat)\(Xmat'*Yvec);

%% instrument construction
xlag            = [X(1:end-1,:)];
xt              = [X(2:end,:)];
n               = nn-h+1;
Rz              = (1-1/(nn^0.95))*eye(m+mf); 
diffx           = xt-xlag; 
z               = zeros(nn,m+mf);
z(1,:)          = diffx(1,:);

for i=2:nn
    z(i,:)=z(i-1,:)*Rz+diffx(i,:);
end
Z               = [zeros(1,m+mf);z(1:n-1,:)];
zz              = [zeros(1,m+mf);z(1:nn-1,:)];
ZK              = zeros(n,m+mf); %% sum over time
for i=1:n
    ZK(i,:)=sum(zz(i:i+h-1,:),1);
end
%construct Zmat
Zmat            = zeros(nn*m,m*(1+mf));
for i = 1:nn
    for j = 1:m
        Zmat((i-1)*m+j,1+(j-1)*(1+mf))  = [Z(i,j)];        
    end       
end
theta_IV        = (Zmat'*Xmat)\(Zmat'*Yvec);

%autoregressive residual estimation 

u =zeros(nn,m+mf);
for i=1:m+mf
%     rn=regress(xt(:,i),[ones(nn,1) xlag(:,i)]);
%     u(:,i) = xt(:,i)-[ones(nn,1) xlag(:,i)]*rn; %
    rn=regress(xt(:,i),[ xlag(:,i)]);
    u(:,i) = xt(:,i)-[ xlag(:,i)]*rn; %
end
%u=xt-xlag*rn;

%covariance matrix estimation (autoregression)
covu=zeros(m+mf,m+mf);
for t=1:nn
    covu=covu+u(t,:)'*u(t,:);
end
covu=covu/nn;

mm=floor(nn^(1/3)); 
uu=zeros(m+mf,m+mf);
for h=1:mm
    a=zeros(m+mf,m+mf);
    for t=(h+1):nn
        a=a+u(t,:)'*u(t-h,:);
    end
    uu=uu+(1-h/(mm+1))*a;
end 
uu=uu/nn;
Omegauu=covu+uu+uu'; 

%% variance-covariance matrix of W0_mat
W0              = Yvec - Xmat*theta;
W0_mat          = reshape(W0,m,nn);
UU              = zeros(m,m);
ZX              = zeros(m*(1+mf),m*(1+mf));
XZ              = zeros(m*(1+mf),m*(1+mf));
ZZ              = zeros(m*(1+mf),m*(1+mf));
ZS              = zeros(m,m*(1+mf));
XX              = zeros(m*(1+mf),m*(1+mf));
XS              = zeros(m*(1+mf),m*(1+mf));
for ii = 1:nn
    UU              = UU + W0_mat(:,ii)*W0_mat(:,ii)';
end
UU              = UU/nn;
covepshat       = UU;

%% covariance matrix between 'W0' and 'u'
covuhat=zeros(m,m+mf);
for i=1:m+mf
    for j=1:m
        covuhat(j,i)=sum(W0_mat(j,:)*u(:,i));
    end
end
covuhat=covuhat'/nn;

residue = zeros(m,m+mf);
for i = 1:m
    q=zeros(mm,m+mf);
    for h=1:mm
        p=zeros(nn-h,m+mf);
            for t=(h+1):nn
                    p(t-h,:)=u(t,:)*W0_mat(i,t-h)'; %u seems to have a problem
            end
        q(h,:)=(1-h/(1+mm))*sum(p);
    end
    residue(i,:) = sum(q)/nn;
end
Omegaeu=covuhat+residue';  %residue is different
FM=covepshat-Omegaeu'*inv(Omegauu)*Omegaeu;

for ii = 1:nn
    ZX              = ZX + (Zmat((ii-1)*m+1:ii*m,:))'*Xmat((ii-1)*m+1:ii*m,:);
    XZ              = XZ + (Xmat((ii-1)*m+1:ii*m,:))'*Zmat((ii-1)*m+1:ii*m,:);
    ZZ              = ZZ + (Zmat((ii-1)*m+1:ii*m,:))'*UU*(Zmat((ii-1)*m+1:ii*m,:));
    ZS              = ZS + (Zmat((ii-1)*m+1:ii*m,:));
    XS              = XS + (Xmat((ii-1)*m+1:ii*m,:))'*UU*(Xmat((ii-1)*m+1:ii*m,:));
    XX              = XX + (Xmat((ii-1)*m+1:ii*m,:))'*Xmat((ii-1)*m+1:ii*m,:);
end
M0                  = ZZ - ZS'*FM*ZS/nn;
V0                  = inv(ZX)*M0*inv(XZ);
V1                  = inv(XX)*XS*inv(XX);

%% construct the IV-based test
betas               = theta_IV(1:1:end);
M                   = max(betas.^2);
V                   = V0;
R1                  = kron(eye(m),[1]);
V_sr                = chol(R1*V*R1')';
AB                  = randn(m,B);
Psi                 = zeros(B,1);

for ii = 1:B
    temp = V_sr*AB(:,ii);
    Psi(ii,1) = max(temp.^2);
end

q  = quantile(Psi,1-epsilon);

if M>q
    rej = 1;
else
    rej = 0;
end

% J-test
J   = betas'*(inv(V))*betas;
pv  = 1-cdf('chi2',J, m);

%% construct the OLS-based F-test
[rej1,F] = F_test(R(2:end),X(1:end-1,:));

%% construct the OLS-based Wald test
[rej_W, W] = Wald_test(R(2:end,1),X(1:end-1,:));

%% construct the OLS-based Wald test (my own program)
[pv_W0, pv_W1] = Wald_test_0429(R(2:end),X(1:end-1,:));

end

%% KMS's original test
function [Aols,Aivx,Wivx,WivxInd,Q,corrmat]=ivxlh_0404(yt,xt,K)
xlag=xt(1:end-1,:);
xt=xt(2:end,:);
y=yt(2:end,:);

[nn,l]=size(xlag); 
X=[ones(nn,1) xlag];

Wivx=zeros(2,1);
WivxInd=zeros(2,l);

%predictive regression residual estimation 
[Aols,bhat,epshat]=regress(y,X); 

rn=zeros(l,l);
for i=1:l
    rn(i,i)=regress(xt(:,i),xlag(:,i));
end

%autoregressive residual estimation 
u=xt-xlag*rn;

%residuals' correlation matrix
corrmat=corrcoef([epshat u]);

%covariance matrix estimation (predictive regression)
covepshat=epshat'*epshat/nn;
covu=zeros(l,l);
for t=1:nn
    covu=covu+u(t,:)'*u(t,:);
end

%covariance matrix estimation (autoregression)
covu=covu/nn;
covuhat=zeros(1,l);
for i=1:l
    covuhat(1,i)=sum(epshat.*u(:,i));
end

%covariance matrix between 'epshat' and 'u'
covuhat=covuhat'/nn; 

m=floor(nn^(1/3)); 
uu=zeros(l,l);
for h=1:m
    a=zeros(l,l);
    for t=(h+1):nn
        a=a+u(t,:)'*u(t-h,:);
    end
    uu=uu+(1-h/(m+1))*a;
end 
uu=uu/nn;
Omegauu=covu+uu+uu'; 

q=zeros(m,l);
for h=1:m
    p=zeros(nn-h,l);
    for t=(h+1):nn
        p(t-h,:)=u(t,:)*epshat(t-h)';
    end
    q(h,:)=(1-h/(1+m))*sum(p);
end
residue=sum(q)/nn;
Omegaeu=covuhat+residue'; 

%instrument construction
n=nn-K+1;
Rz=(1-1/(nn^0.95))*eye(l); 
diffx=xt-xlag; 
z=zeros(nn,l);
z(1,:)=diffx(1,:);
for i=2:nn
    z(i,:)=z(i-1,:)*Rz+diffx(i,:);
end
Z=[zeros(1,l);z(1:n-1,:)];


zz=[zeros(1,l);z(1:nn-1,:)];
ZK=zeros(n,l);
for i=1:n
    ZK(i,:)=sum(zz(i:i+K-1,:),1);
end

yy=zeros(n,1);
for i=1:n
    yy(i)=sum(y(i:i+K-1));
end 
xK=zeros(n,l);
for i=1:n
    xK(i,:)=sum(xlag(i:i+K-1,:),1);
end 

meanxK=mean(xK);
Yt=yy-mean(yy);
Xt=zeros(n,l);
for i=1:l
    Xt(:,i)=xK(:,i)-meanxK(:,i)*ones(n,1);
end

Aivx=Yt'*Z*pinv(Xt'*Z);
meanzK=mean(ZK);

FM=covepshat-Omegaeu'*Omegauu^(-1)*Omegaeu;
M=ZK'*ZK*covepshat-n*meanzK'*meanzK*FM;

H=eye(l);
Q=H*pinv(Z'*Xt)*M*pinv(Xt'*Z)*H';
Wivx(1,1)=(H*Aivx')'*pinv(Q)*(H*Aivx');   
Wivx(2,1)= 1-cdf('chi2',Wivx(1,1), l);
    
WivxInd(1,:)=(Aivx./((diag(Q)).^(1/2))').^2;
WivxInd(2,:)=1-cdf('chi2',WivxInd(1,:), 1);
    
end


function [pv,F] = F_test(Y,XX);
% this function is to test if all coefficients (except for the constant)
% are zero

[n p]   = size(XX);
con     = ones(n,1);
XX      = [con XX];
beta    = (inv(XX'*XX)*(XX'*Y));
SSE_full= sum((Y-XX*beta).^2);
df_full = n-p-1;

beta1   = (inv(con'*con)*(con'*Y));
SSE_red = sum((Y-con*beta1).^2);

F       = ((SSE_red-SSE_full)/p)/(SSE_full/df_full);
pv      = 1-fcdf(F,p,df_full);
end

function [Aols,Aivx,Wivx,WivxInd,Q,corrmat]=ivxlh_modified_0419(yt,xt,K,printres)
xlag=xt(1:end-1,:);
xt=xt(2:end,:);
y=yt(2:end,:);

[nn,l]=size(xlag); 
X=[ones(nn,1) xlag];

Wivx=zeros(2,1);
WivxInd=zeros(2,l);

%predictive regression residual estimation 
[Aols,bhat,epshat]=regress(y,X); 

[~,~,epshat] = regress(y,ones(nn,1)); %%major change

rn=zeros(l,l);
for i=1:l
    rn(i,i)=regress(xt(:,i),xlag(:,i));
end

%autoregressive residual estimation 
u=xt-xlag*rn;

%residuals' correlation matrix
corrmat=corrcoef([epshat u]);

%covariance matrix estimation (predictive regression)
covepshat=epshat'*epshat/nn;
covu=zeros(l,l);
for t=1:nn
    covu=covu+u(t,:)'*u(t,:);
end

%covariance matrix estimation (autoregression)
covu=covu/nn;
covuhat=zeros(1,l);
for i=1:l
    covuhat(1,i)=sum(epshat.*u(:,i));
end

%covariance matrix between 'epshat' and 'u'
covuhat=covuhat'/nn; 

m=floor(nn^(1/3)); 
uu=zeros(l,l);
for h=1:m
    a=zeros(l,l);
    for t=(h+1):nn
        a=a+u(t,:)'*u(t-h,:);
    end
    uu=uu+(1-h/(m+1))*a;
end 
uu=uu/nn;
Omegauu=covu+uu+uu'; 

q=zeros(m,l);
for h=1:m
    p=zeros(nn-h,l);
    for t=(h+1):nn
        p(t-h,:)=u(t,:)*epshat(t-h)';
    end
    q(h,:)=(1-h/(1+m))*sum(p);
end
residue=sum(q)/nn;
Omegaeu=covuhat+residue'; 

%instrument construction
n=nn-K+1;
Rz=(1-1/(nn^0.95))*eye(l); 
diffx=xt-xlag; 
z=zeros(nn,l);
z(1,:)=diffx(1,:);
for i=2:nn
    z(i,:)=z(i-1,:)*Rz+diffx(i,:);
end
Z=[zeros(1,l);z(1:n-1,:)];


zz=[zeros(1,l);z(1:nn-1,:)];
ZK=zeros(n,l);
for i=1:n
    ZK(i,:)=sum(zz(i:i+K-1,:),1);
end

yy=zeros(n,1);
for i=1:n
    yy(i)=sum(y(i:i+K-1));
end 
xK=zeros(n,l);
for i=1:n
    xK(i,:)=sum(xlag(i:i+K-1,:),1);
end 

meanxK=mean(xK);
Yt=yy-mean(yy);
Xt=zeros(n,l);
for i=1:l
    Xt(:,i)=xK(:,i)-meanxK(:,i)*ones(n,1);
end

Aivx=Yt'*Z*pinv(Xt'*Z);
meanzK=mean(ZK);

FM=covepshat-Omegaeu'*Omegauu^(-1)*Omegaeu;
M=ZK'*ZK*covepshat-n*meanzK'*meanzK*FM;

H=eye(l);
Q=H*pinv(Z'*Xt)*M*pinv(Xt'*Z)*H';
Wivx(1,1)=(H*Aivx')'*pinv(Q)*(H*Aivx');   
Wivx(2,1)= 1-cdf('chi2',Wivx(1,1), l);
    
WivxInd(1,:)=(Aivx./((diag(Q)).^(1/2))').^2;
WivxInd(2,:)=1-cdf('chi2',WivxInd(1,:), 1);
    
end

function [pv, W] = Wald_test(Y,X)
[n m_x] = size(X);
X = [ones(n,1) X]; %add a constant
b = (inv(X'*X)*(X'*Y));
sig2 = (sum((Y-X*b).^2))/(n-1-(m_x+1));
V = (inv(X'*X)*sig2);
W = (b(2:end,1)'*inv(V(2:end,2:end))*b(2:end,1));
pv = 1-cdf('chi2',W, m_x);
end

function [pv0, pv1] = Wald_test_0429(Y,X)
[n m_x] = size(X);
X = [ones(n,1) X]; %add a constant
b = (inv(X'*X)*(X'*Y));
res = (Y-X*b).^2;
temp=zeros(m_x+1,m_x+1);
for i = 1:n
    temp=temp+ (X(i,:)'*X(i,:))*res(i);
end
Omega = temp/n;
EXX   = inv(X'*X/n);
V0 = (EXX*Omega*EXX)/n;
W0 = (b(2:end,1)'*inv(V0(2:end,2:end))*b(2:end,1));
pv0 = 1-cdf('chi2',W0, m_x);
V1 = (n/(n-(m_x+1)))*V0;
W1 = (b(2:end,1)'*inv(V1(2:end,2:end))*b(2:end,1));
pv1 = 1-cdf('chi2',W1, m_x);
end

% This function is used to compute OLS-F test (from Cochrane's website)
function [bv,sebv,R2v,R2vadj,v,F] = olsgmm(lhv,rhv,lags,weight);

if size(rhv,1) ~= size(lhv,1);
   disp('olsgmm: left and right sides must have same number of rows. Current rows are');
   size(lhv)
   size(rhv)
end;

T = size(lhv,1);
N = size(lhv,2);
K = size(rhv,2);
sebv = zeros(K,N);
Exxprim = inv((rhv'*rhv)/T);
bv = rhv\lhv;

% if weight, lags are is scalar expand so all have the same value
if (size(weight,1) == 1) & (size(lhv,2) > 1); 
    weight = weight*ones(size(lhv,2),1); 
end; 
if (size(lags,1) == 1) & (size(lhv,2) > 1); 
    lags = lags*ones(size(lhv,2),1); 
end; 


if weight == -1;  % skip ses if you don't want them.  returns something so won't get error message
    sebv=NaN;
    R2v=NaN;
    R2vadj=NaN;
    v=NaN;
    F=NaN;
else; 
    errv = lhv-rhv*bv;
    s2 = mean(errv.^2);
    vary = lhv - ones(T,1)*mean(lhv);
    vary = mean(vary.^2);

    R2v = (1-s2./vary)';
    R2vadj= (1 - (s2./vary)*(T-1)/(T-K))';
    
    %compute standard errors
    for indx = 1:N;
        err=errv(:,indx);
        if (weight(indx) == 0)|(weight(indx) == 1)
            inner = (rhv.*(err*ones(1,K)))'*(rhv.*(err*ones(1,K)))/T;
            
        	for jindx = (1:lags(indx));
                inneradd = (rhv(1:T-jindx,:).*(err(1:T-jindx)*ones(1,K)))'...
        	              *(rhv(1+jindx:T,:).*(err(1+jindx:T)*ones(1,K)))/T;
        	    inner = inner + (1-weight(indx)*jindx/(lags(indx)+1))*(inneradd + inneradd');
            end;
        elseif weight(indx) == 2; 
            inner = rhv'*rhv/T; 
            for jindx = 1:lags(indx); 
                inneradd = rhv(1:T-jindx,:)'*rhv(1+jindx:T,:)/T;
                inner = inner + (1-jindx/lags(indx))*(inneradd+inneradd'); 
            end; 
            inner = inner*std(err)^2;
        end; 
        varb = 1/T*Exxprim*inner*Exxprim;
        
        % F test for all coeffs (except constant) zero -- actually chi2 test
        if rhv(:,1) == ones(size(rhv,1),1); 
            chi2val = bv(2:end,indx)'*inv(varb(2:end,2:end))*bv(2:end,indx);
            dof = size(bv(2:end,1),1); 
            pval = 1-cdf('chi2',chi2val, dof); 
            F(indx,1:3) = [chi2val dof pval]; 
        else; 
            chi2val = bv(:,indx)'*inv(varb)*bv(:,indx);
            dof = size(bv(:,1),1); 
            pval = 1-cdf('chi2',chi2val, dof); 
            F(indx,1:3) = [chi2val dof pval]; 
        end; 
            
        if indx == 1; 
           v = varb;
        else;
           v = [v; varb ];
        end;
        
       seb = diag(varb);
       seb = sign(seb).*(abs(seb).^0.5);
       sebv(:,indx) = seb;
    end;
end; % ends if w > -1;     

end

```

## Monte Carlo Simulation Studies:

Here we present the R coding and some preliminary results based on a Monte Carlo Simulation study of the finite-sample properties of the inferential procedure that is under development for this research project.  

Consider the following Toeplitz structure of the covariance matrix. 



```R



```

## Main References:

- Katsouris, C. (2022c). "Estimation and Inference in Quantile Predictive Regression Systems" (Chapter 4, PhD thesis, School of Economic, Social and Political Sciences, University of Southampton.
- Katsouris, C. (2022b). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregression and Predictive Regression Models". University of Southampton, Department of Economics, Working paper.  
- Katsouris, C. (2022a). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregressive Time Series". arXiv preprint [arXiv:2204.02073](https://arxiv.org/abs/2204.02073).
- Katsouris, C. (2021e). "Bootstrapping Nonstationary Autoregressive Processes in Predictive Regression". University of Southampton, Working paper.   
- Katsouris, C. (2021d). "Testing for Structural Breaks in Predictive Regression Models". University of Southampton, Department of Economics, Working paper.  
- Ke-Li Xu & Junjie Guo (2022). "A New Test for Multiple Predictive Regression". Journal of Financial Econometrics.
- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). "Robust econometric inference for stock return predictability". The Review of Financial Studies, 28(5), 1506-1553.
- Koo, B., Anderson, H. M., Seo, M. H., & Yao, W. (2020). "High-dimensional predictive regression in the presence of cointegration". Journal of Econometrics, 219(2), 456-477.
- Lee, J. H. (2016). "Predictive quantile regression with persistent covariates: IVX-QR approach". Journal of Econometrics, 192(1), 105-118.
- Lee, J. H., Shi, Z., & Gao, Z. (2022). "On LASSO for predictive regression". Journal of Econometrics, 229(2), 322-349.
- Fan, R., & Lee, J. H. (2019). "Predictive quantile regressions under persistence and conditional heteroskedasticity". Journal of Econometrics, 213(1), 261-280.
- Chen, W. W., Deo, R. S., & Yi, Y. (2013). "Uniform inference in predictive regression models". Journal of Business & Economic Statistics, 31(4), 525-533.
- Magdalinos, T., and Petrova, K. (2022). "Uniform and distribution-free inference with general autoregressive processes". University of Southampton, Department of Economics, Working paper.
- Magdalinos, T. (2016). "Least squares and IVX limit theory in systems of predictive regressions with GARCH innovations". Econometric Theory, 1-38.
- Magdalinos, T., and Phillips, P. C. B. (2009). "Limit theory for cointegrated systems with moderately integrated and moderately explosive regressors". Econometric Theory, 25(2), 482-526.
- Phillips, P. C. B., and Magdalinos, T. (2009). "Econometric inference in the vicinity of unity". Singapore Management University, CoFie Working paper, 7.
- Phillips, P. C. B., and Magdalinos, T. (2008). "Limit theory for explosively cointegrated systems". Econometric Theory, 24(4), 865-887.
- Phillips, P. C. B., and Magdalinos, T. (2007). "Limit theory for moderate deviations from a unit root". Journal of Econometrics, 136(1), 115-130.
- Wagner, M., Grabarczyk, P., & Hong, S. H. (2020). "Fully modified OLS estimation and inference for seemingly unrelated cointegrating polynomial regressions and the environmental Kuznets curve for carbon dioxide emissions". Journal of Econometrics, 214(1), 216-255.
- Yousuf, K., & Ng, S. (2021). "Boosting high dimensional predictive regressions with time varying parameters". Journal of Econometrics, 224(1), 60-87.
- Zhu, F., Cai, Z., & Peng, L. (2014). "Predictive regressions for macroeconomic data". The Annals of Applied Statistics, 8(1), 577-594.

## Literature on high-dimensional IV regression:

- Krampe, J., Paparoditis, E., & Trenkler, C. (2022). Structural inference in sparse high-dimensional vector autoregressions. Journal of Econometrics.
- Belloni, A., Hansen, C., & Newey, W. (2017). Simultaneous confidence intervals for high-dimensional linear models with many endogenous variables. arXiv preprint arXiv:1712.08102.
- Gold, D., Lederer, J., & Tao, J. (2020). "Inference for high-dimensional instrumental variables regression". Journal of Econometrics, 217(1), 79-111.
- Ning, Y., & Liu, H. (2017). A general theory of hypothesis tests and confidence regions for sparse high dimensional models. The Annals of Statistics, 45(1), 158-195.
- Van de Geer, S., Bühlmann, P., Ritov, Y. A., & Dezeure, R. (2014). On asymptotically optimal confidence regions and tests for high-dimensional models. The Annals of Statistics, 42(3), 1166-1202.

# Acknowledgments

The author greatfully acknowledges financial support from the [Department of Economics](http://business-school.exeter.ac.uk/about/departments/economics/) of the [Faculty of Environment, Science and Economy](https://www.exeter.ac.uk/departments/ese/) at the University of Exeter, United Kingdom. 

Christis G. Katsouris is a Lecturer in Economics at the [University of Exeter Business School](http://business-school.exeter.ac.uk/). He is also a member of the [Time Series and Machine Learning Group](https://www.personal.soton.ac.uk/cz1y20/Reading_Group/mlts-group-2022.html) at the [School of Mathematical Sciences](https://www.southampton.ac.uk/about/faculties-schools-departments/school-of-mathematical-sciences) (Statistics Division) of the University of Southampton. 

# Declarations

The author (Christis G. Katsouris) declares no conflicts of interest.

Notice that the academic research presented here is considered to be as open access to the academic and non-academic community. Therefore, we would appreciate it if appropriate acknolwedgement is given to statistical methodologies and econometric procedures developed by academic researchers and made available to the wider applied data scientist community.   

# Historical Background

#### Harald Cramér 

Harald Cramér  (25 September 1893 – 5 October 1985) was a Swedish mathematician, actuary, and statistician, specializing in mathematical statistics and probabilistic number theory. John Kingman described him as "one of the giants of statistical theory". A large portion of Cramér's work concerned the field of actuarial science and insurance mathematics. In 1929, Cramér was appointed to a newly created chair in Stockholm University, becoming the first Swedish professor of Actuarial Mathematics and Mathematical Statistics. Cramér retained this position up until 1958. During his tenure at Stockholm University, Cramér was a PhD advisor for 10 students, most notably Herman Wold and Kai Lai Chung. In 1950 he was elected as a Fellow of the American Statistical Association. Starting in 1950, Cramér took on the additional responsibility of becoming the President of Stockholm University. In 1958, he was also appointed to be Chancellor of the entire Swedish university system (Source: Wikepedia). 

<p align="center">
  
<img src="https://github.com/christiskatsouris/Robust_high_dimensional_IV_Cointegration/blob/main/data/BTC.jpg" width="560"/>

</p>  

> A two-year explosive bubble has just burst on the 13th of November 2022!
