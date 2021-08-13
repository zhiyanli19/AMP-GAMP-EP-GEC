function [MSE_error] = EPnew(obj, Input)
eps = 5e-7;
IterNum=Input.IterNum;
%mes = Input.mes;

N = Input.N;
M = Input.M;
nuw = Input.nuw;
H = obj.H;
xo = obj.xo;
x = obj.x;
y = obj.y1;



mes=Input.mes;
IterNum=Input.IterNum;
m0=zeros(N,1);
v0=ones(N,1);
v0_inv=1./v0;
v1=ones(M,1);
v1_inv=1./v1;
m1=zeros(M,1);
v1_inv_old=v1_inv;
m1_old=m1;
m0_sub = zeros(N, 1);
m1_sub = ones(N, 1);
v0_inv_old=v0_inv;
m0_old=m0;
m1_sub_old = m1_sub;
m0_sub_old = m0_sub;
MSE_error=zeros(IterNum,1);
MSE_old=1;
for ii=1:IterNum
  
    [hatx_plus,vx_plus]=estimator_x(xo,m0,v0);

    vx_plus=max(vx_plus,5e-13);
    v1_inv=(v0-vx_plus)./vx_plus./v0;
    m1_sub=(hatx_plus.*v0-m0.*vx_plus)./vx_plus./v0;
    
    negldx=v1_inv<eps;
    v1_inv(negldx)=v1_inv_old(negldx);
    m1_sub(negldx)=m1_sub_old(negldx);
    
    % damping
    [v1_inv,v1_inv_old]=damping(v1_inv,v1_inv_old, mes);
    [m1_sub,m1_sub_old]=damping(m1_sub,m1_sub_old, mes);
      
    Qx_plus=(H'*diag(1 / nuw)*H+diag(v1_inv))^(-1);
    mx_plus=Qx_plus*(H'*diag(1 / nuw)*y+m1_sub);
    vz_plus=real(diag(Qx_plus));
    
    MSE=norm(mx_plus-x).^2/N;
    if sum(isnan(vz_plus))>0 || MSE>MSE_old
        MSE_error(ii:IterNum)=MSE_old;
        break;
    else
        MSE_old=MSE;
    end 
    MSE_error(ii,1)=MSE;
    v1 = 1./v1_inv;
    m1 = m1_sub .* v1;
    vz_plus=max(vz_plus,5e-13);
    v0_inv=(v1-vz_plus)./v1./vz_plus;
    m0_sub=(mx_plus.*v1-m1.*vz_plus)./v1./vz_plus;
    %% ÌÞ³ý¸ºÔªËØ 
    negldx=v0_inv<0;
    v0_inv(negldx)=v0_inv_old(negldx);
    m0_sub(negldx)=m0_sub_old(negldx);
     
     %% damping
    [v0_inv,v0_inv_old]=damping(v0_inv,v0_inv_old,mes);
    [m0_sub,m0_sub_old]=damping(m0_sub,m0_sub_old,mes);
    v0 = 1 ./ v0_inv;  
    m0 = m0_sub .* v0;
end
end



function [umean,uvar]=estimator_x(xo,v,wvar)
logpxr = bsxfun(@times, -1./wvar, abs(bsxfun(@minus, v, xo)).^2);
% logpxr = bsxfun(@minus, logpxr, max(logpxr) );  % for stability          
pxr = exp(logpxr);
pxr = bsxfun(@rdivide, pxr, sum(pxr,2) );
umean = sum(bsxfun(@times, pxr, xo), 2); 
uvar = sum(pxr .* abs(bsxfun(@minus, umean, xo)).^2, 2);
end

function [x,x_old]=damping(x, x_old, mes)
x=mes*x+(1-mes)*x_old;
x_old=x;
end
