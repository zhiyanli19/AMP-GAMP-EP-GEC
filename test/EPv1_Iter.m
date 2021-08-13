function [ EP_MSE] = EPv1_Iter(obj, Input)
    eps = 5e-13;  
    N = Input.N;
    M = Input.M;
    H = obj.H;
    xo = obj.xo;
    x= obj.x;
    y = obj.y;
    nuw = Input.nuw;
    %% ³õÊ¼»¯
    mx = zeros( N , 1 );
    mxold=mx;
    vx = ones( N , 1 );
    vxinv=1./vx;
    vxinvold=vxinv;
    L = Input.iter_num;
    EP_MSE = zeros( 1 , L );
    MSE_old=1;
   %% µü´ú
   for ii=1:L
       Vy_hat=(diag(vxinv) + H'*(1./nuw) * H )^(-1);
       my_hat= Vy_hat * ( mx + H'*(1./nuw)*y  );
       
       vy_hat=real( diag(Vy_hat) );
       %vy=1./ ( 1./vy_hat-1./vx );
        vy=max(vy_hat./(1-vy_hat.*vxinv),eps);
        my=vy.*( my_hat./vy_hat - mx) ;   
        [mx_hat, vx_hat] = Estimator_1( xo, my, vy );
        vx_hat = max( vx_hat , eps);
        vxinv=(vy-vx_hat)./vx_hat./vy;
        %vx=1./(1./vx_hat-1./vy);
        mx=(mx_hat.*vy-my.*vx_hat)./vx_hat./vy;
        
        MSE = norm (my_hat - x ).^2 / N;
        if sum(isnan(vy_hat))>0 || MSE>MSE_old
            EP_MSE(ii-0:L)=MSE_old;
            break;
        else
            MSE_old=MSE;
        end
        EP_MSE(1,ii)=MSE;  
        p = vxinv < 0;
        mx(p) = mxold(p);
        vxinv(p) = vxinvold(p);
        [mx,mxold]= damping(mx, mxold, 0.9);
        [vxinv,vxinvold] = damping(vxinv, vxinvold, 0.9);
   end
end
function [umean,uvar]=Estimator_1(xo,v,wvar)
logpxr = bsxfun(@times, -1./wvar, abs(bsxfun(@minus, v, xo)).^2);
logpxr = bsxfun(@minus, logpxr, max(logpxr) );            
pxr = exp(logpxr);
pxr = bsxfun(@rdivide, pxr, sum(pxr,2) );
umean = sum(bsxfun(@times, pxr, xo), 2); 
uvar = sum(pxr .* abs(bsxfun(@minus, umean, xo)).^2, 2);
end
function [x,x_old]=damping(x, x_old, mes)
x=mes*x+(1-mes)*x_old;
x_old=x;
end