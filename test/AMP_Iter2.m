function [MSE_error] = AMP_Iter(obj, Input)
    eps = 1e-30;
    %mes = Input.mes;
    N = Input.N;
    M = Input.M;
    nuw = Input.nuw;
    iter_num = Input.IterNum;
    H = obj.H;
    xo = obj.xo;
    x = obj.x;
    y1 = obj.y1;

    %% array initialization
    hat_x = zeros(N, 1);
    hat_v = ones(N, 1);
    sqrH = abs(H).^2;
    sqrHt = sqrH';
    Ht = H';
    %V = zeros(M, 1);
    Z = zeros(M, 1);
    Z_old=Z;
    sigma = nuw * ones(M, 1);  
    MSE_error = zeros(1, iter_num);    
    hat_x_old = hat_x;
    hat_v_old = hat_v;
    MSE_old=1;
    for ii = 1: iter_num
        %% Output Nodes
        V = sqrH * hat_v;
        %V = max ( V ,eps );
        Z = H * hat_x -  V.*(y1-Z)./(sigma + V );
        [Z, Z_old] = damping(Z, Z_old, 0.8);
        %[V, V_old] = damping(V, V_old, 0.8);
        %% Input Nodes
        vitiled = ( sqrHt * ( (sigma+V).^(-1)) ).^(-1);
        %vitiled = max ( vitiled,eps );
        mitiled = hat_x + vitiled .* (Ht * ((y1-Z)./(sigma+V)));

       [hat_x, hat_v] = Estimator_1(xo, mitiled, vitiled);

        MSE = norm (hat_x - x ).^2 / N;
        %MSE = mean(abs(hat_x - x).^2);
        if sum(isnan(hat_v))>0 || MSE>MSE_old
            MSE_error(ii-1:iter_num)=MSE_old;
            break;
        else
            MSE_old=MSE;
        end
        MSE_error(1, ii) = MSE;
        %hat_v = max(hat_v, eps);
        %x_hat = hat_x(1:N/2) + 1j*hat_x(K/2+1:K);
        [hat_x, hat_x_old] = damping(hat_x, hat_x_old, 0.8);
        [hat_v, hat_v_old] = damping(hat_v, hat_v_old, 0.8);

        %MSE_error(1, ii) = mean(abs(hat_x - obj.x).^2);
    end
    %save('hat_x.mat', 'hat_x');
    %[BER, SER] = BER_Calculation(hat_x, obj, Input);
end

function [x,x_old]=damping(x, x_old, mes)
x=mes*x+(1-mes)*x_old;
x_old=x;
end

   
function [umean,uvar]=Estimator_1(xo,v,wvar)
logpxr = bsxfun(@times, -1./wvar, abs(bsxfun(@minus, v, xo)).^2);
logpxr = bsxfun(@minus, logpxr, max(logpxr) );            
pxr = exp(logpxr);
pxr = bsxfun(@rdivide, pxr, sum(pxr,2) );
umean = sum(bsxfun(@times, pxr, xo), 2); 
uvar = sum(pxr .* abs(bsxfun(@minus, umean, xo)).^2, 2);
end

function [m, v] = estimator_1(check_m, check_v)
    m = zeros(size(check_m));
    v = zeros(size(check_m));
    syms x_sys v_sys m_sys;
    tmp1 = int(x_sys * exp(-(x_sys - m_sys)^2 / (2 * v_sys)), x_sys, 0, 1);
    tmp2 = int(exp(-(x_sys - m_sys)^2 / (2 * v_sys)), x_sys, 0, 1);
    %eps = 1e-50;
    %eta = tmp1 / (tmp2 + eps);
    eta = tmp1 / tmp2;
    tmp3 = int(x_sys^2 * exp(-(x_sys - m_sys)^2 / (2 * v_sys)), x_sys, 0, 1);
    %theta = tmp3 / (tmp2 + eps) - eta^2;
    theta = tmp3 / tmp2 - eta^2;

    for i = 1: size(check_m, 1)
        m(i) = double(subs(eta, [m_sys, v_sys], [check_m(i), check_v(i)]));
        v(i) = double(subs(theta, [m_sys, v_sys], [check_m(i), check_v(i)]));
    end
end
