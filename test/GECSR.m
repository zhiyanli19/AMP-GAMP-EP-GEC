function MSE_error=GECSR(obj, Input)

H=obj.H;
xo=obj.xo;
x=obj.x;

[M,N]=size(H);
mes=Input.mes;
IterNum=Input.IterNum;

v1_plus=ones(M,1);
v1_plus_inv=1./v1_plus;
m1_plus=zeros(M,1);
v1_plus_inv_old=v1_plus_inv;
m1_plus_old=m1_plus;

v0_plus=2*ones(N,1);
v0_plus_inv=1./v0_plus;
m0_plus=zeros(N,1);
v0_plus_inv_old=v0_plus_inv;
m0_plus_old=m0_plus;

MSE_error=zeros(IterNum,1);
MSE_old=1;
for ii=1:IterNum
    
    
    [hatz_sub,vz_sub]=estimator_z(m1_plus./v1_plus_inv,1./v1_plus_inv,obj,Input);
    v1_sub=max(vz_sub./(1-vz_sub.*v1_plus_inv),eps);          %避免出现负数以及0
    m1_sub=v1_sub.*(hatz_sub./vz_sub-m1_plus);
        
    Qx_sub=(H'*diag(1./v1_sub)*H+diag(v0_plus_inv))^(-1);
    hatx_sub=Qx_sub*(H'*diag(1./v1_sub)*m1_sub+m0_plus);                     %这里的r0_plus
    vx_sub=real(diag(Qx_sub));
    v0_sub=max(vx_sub./(1-vx_sub.*v0_plus_inv),eps);
    m0_sub=v0_sub.*(hatx_sub./vx_sub-m0_plus);
    
    %% Forward passing
    [hatx_plus,vx_plus]=estimator_x(xo,m0_sub,v0_sub);
    MSE=norm(hatx_plus-x).^2/N;
    if sum(isnan(vx_plus))>0 || MSE>MSE_old
        MSE_error(ii-1:IterNum)=MSE_old;
        break;
    else
        MSE_old=MSE;
    end 
    MSE_error(ii,1)=MSE;
    vx_plus=max(vx_plus,5e-13);
    v0_plus_inv=(v0_sub-vx_plus)./vx_plus./v0_sub;
    m0_plus=(hatx_plus.*v0_sub-m0_sub.*vx_plus)./vx_plus./v0_sub;
    
    negldx=v0_plus_inv<eps;
    v0_plus_inv(negldx)=v0_plus_inv_old(negldx);
    m0_plus(negldx)=m0_plus_old(negldx);
    
    % damping
    [v0_plus_inv,v0_plus_inv_old]=damping(v0_plus_inv,v0_plus_inv_old, mes);
    [m0_plus,m0_plus_old]=damping(m0_plus,m0_plus_old, mes);
      
    Qx_plus=(H'*diag(1./v1_sub)*H+diag(v0_plus_inv))^(-1);
    mx_plus=Qx_plus*(H'*diag(1./v1_sub)*m1_sub+m0_plus);
    hatz_plus=H*mx_plus;
    vz_plus=real(diag(H*Qx_plus*H'));
          
    vz_plus=max(vz_plus,5e-13);
    v1_plus_inv=(v1_sub-vz_plus)./v1_sub./vz_plus;
    m1_plus=(hatz_plus.*v1_sub-m1_sub.*vz_plus)./v1_sub./vz_plus;
    %% 剔除负元素 
    negldx=v1_plus_inv<0;
    v1_plus_inv(negldx)=v1_plus_inv_old(negldx);
    m1_plus(negldx)=m1_plus_old(negldx);
     
     %% damping
    [v1_plus_inv,v1_plus_inv_old]=damping(v1_plus_inv,v1_plus_inv_old,mes);
    [m1_plus,m1_plus_old]=damping(m1_plus,m1_plus_old,mes);
      
    
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


function [hatz,hatv]=estimator_z(m,v,obj,Input)
nuw=Input.nuw;
ADC_switch=Input.ADC_switch;
y2=obj.y2;
M=Input.M;
DeltaTh=obj.DeltaTh;
quan_step=obj.quan_step;

if ADC_switch==0
    hatv=1./(1./v+1/nuw);
    hatz=hatv.*(y2/nuw+m./v);
else
    y2=[real(y2);imag(y2)];
    m=[real(m);imag(m)];
    v=[real(v);real(v)];
    y_up=y2+quan_step/2;
    y_low=y2-quan_step/2;
    [pos1,~]=find(y2>max(DeltaTh));   
    [pos2,~]=find(y2<-max(DeltaTh));
    y_up(pos1)=1e5;
    y_low(pos2)=-1e5;

    
    eta1=(y_up-m)./sqrt((nuw+v)/2);
    eta2=(y_low-m)./sqrt((nuw+v)/2);
      
    tem1=normpdf(eta1)-normpdf(eta2);
    tem2=normcdf(eta1)-normcdf(eta2);
    tem3=eta1.*normpdf(eta1)-eta2.*normpdf(eta2);
    
%     tem1(pos1)=-normpdf(eta2(pos1));
%     tem1(pos2)=normpdf(eta1(pos2));
%     tem2(pos1)=1-normcdf(eta2(pos1));
%     tem2(pos2)=normcdf(eta1(pos2));
%     tem3(pos1)=-eta2(pos1).*normpdf(eta2(pos1));
%     tem3(pos2)=eta1(pos2).*normpdf(eta1(pos2));
    
    z_tem=m-v./(sqrt(2*(nuw+v))).*(tem1./tem2);
    v_tem=v/2-((v.^2)./(2*(nuw+v))).*(tem3./tem2+(tem1./tem2).^2);
    
    hatz=z_tem(1:M)+1j*z_tem(M+1:2*M);
    hatv=max(v_tem(1:M)+v_tem(M+1:2*M),eps);         
end
end