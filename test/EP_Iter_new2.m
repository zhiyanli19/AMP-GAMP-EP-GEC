function [MSE_error] = EP_Iter_new2(obj, Input)
    eps = 5e-7;
    IterNum=Input.IterNum;
    %mes = Input.mes;

    N = Input.N;
    M = Input.M;
    nuw = Input.nuw;
    H = obj.H;
    xo = obj.xo;
    x = obj.x;
    y1 = obj.y1;

    %% array initialization
    r = zeros(N, 1);
    A = ones(N, 1);
    Ht = H';



    MSE_error=zeros(IterNum,1);
    MSE_old=1;


    for ii = 1: IterNum

        V = inv(1./nuw * H' * H + diag(A));
        m = V * (1 ./ nuw * H' * y1 + r);

        V_i = real(diag(V));
        h_i = V_i ./ (1 - V_i .* A);
        h_i = max(h_i, eps);
        t_i = h_i .* (m ./ V_i - r);
        [m_pi, V_pi] = estimator_x(xo, t_i, h_i);
        V_pi = max(V_pi, eps);
        A_new = (h_i - V_pi) ./ V_pi ./ h_i;
        r_new = (m_pi .* h_i - t_i .* V_pi) ./ V_pi ./ h_i;
        
       MSE=norm(m-x).^2/N;
       if sum(isnan(V_i))>0 || MSE>MSE_old
           MSE_error(ii-1:IterNum)=MSE_old;
           break;
       else
           MSE_old=MSE;
       end 
       MSE_error(ii,1)=MSE;
        
        position = A_new < 0;
        A_new(position) = A(position);
        r_new(position) = r(position);
        %A = A_new;
        %r = r_new;
        %if sum(isnan(hat_v)) > 0 
            %hat_x = hat_x_old;
            %hat_v = hat_v_old;
            %break;
        %end

        [r_new, r] = damping(r_new, r, 0.8);
        [A_new, A] = damping(A_new, A, 0.8);

    end

    %[SER] = BER_Calculation(m, obj, Input);
    %MSE_error = mean(abs(m - obj.x).^2);
end

function [x, x_old] = damping(x, x_old, mes)
    x = mes * x + (1 - mes) * x_old;
    x_old = x;
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

function [m, v] = Estimator_1(xo, check_m, check_v)
    eps = 1e-30;
    log_posterior = bsxfun(@times, -1 ./ check_v, abs(bsxfun(@minus, xo, check_m).^2));
    log_posterior = bsxfun(@minus, log_posterior, max(log_posterior));  %防止溢出

    posterior = exp(log_posterior); 
    posterior = bsxfun(@rdivide, posterior, sum(posterior, 2));   %得到标准PDF
    m = sum(bsxfun(@times, posterior, xo), 2);                   %计算PDF的均�?
    v = sum(posterior .* abs(bsxfun(@minus, m, xo).^2), 2);       %计算PDF的方�?
end
function [umean,uvar]=estimator_x(xo,v,wvar)
logpxr = bsxfun(@times, -1./wvar, abs(bsxfun(@minus, v, xo)).^2);
% logpxr = bsxfun(@minus, logpxr, max(logpxr) );  % for stability          
pxr = exp(logpxr);
pxr = bsxfun(@rdivide, pxr, sum(pxr,2) );
umean = sum(bsxfun(@times, pxr, xo), 2); 
uvar = sum(pxr .* abs(bsxfun(@minus, umean, xo)).^2, 2);
end

function [hat_s, hat_tau] = estimator_2(Z, V, obj, Input)
    %eps = 1e-30;
    nuw = Input.nuw;
    y = obj.y;
    M = Input.M;
    quan_step = obj.quan_step;
    Quan_bound = (2^(Input.bit - 1) - 1) * quan_step;

    y_ = [real(y); imag(y)];
    z_ = [real(Z); imag(Z)];
    v_ = [real(V); real(V)];

    y_up = y_ + quan_step / 2;
    y_low = y_ - quan_step / 2;
    [pos1, ~] = find(y_up > Quan_bound);
    [pos2, ~] = find(y_low < -Quan_bound);
    threhold = 4 * Quan_bound;
    %threhold = 1e5;
    y_up(pos1) = threhold;
    y_low(pos2) = -threhold;

    eta1 = (sign(y_) .* z_ - min(abs(y_up), abs(y_low))) ./ sqrt((nuw + v_) / 2);
    eta2 = (sign(y_) .* z_ - max(abs(y_up), abs(y_low))) ./ sqrt((nuw + v_) / 2);
    tem1 = normpdf(eta1) - normpdf(eta2);
    tem2 = normcdf(eta1) - normcdf(eta2);
    tem3 = eta1 .* normpdf(eta1) - eta2 .* normpdf(eta2);
    %pos = eta2 < -100; 
    %tem1(pos) = normpdf(eta1(pos));
    %tem2(pos) = normcdf(eta1(pos));
    %tem3(pos) = eta1(pos) .* normpdf(eta1(pos));

    z_tem = z_ + (sign(y_) .* v_ ./ sqrt(2 * (nuw + v_))) .* (tem1 ./ tem2);
    %z_tem = z_ + (sign(y_) .* v_ ./ sqrt(2 * (nuw + v_))) .* (tem1 ./ tem2);
    v_tem = v_ / 2 - (v_.^2 ./ (2 * (nuw + v_))) .* (tem3 ./ tem2 + (tem1 ./ tem2).^2);
    %v_tem = v_ / 2 - ((v_.^2) ./ (2 * (nuw + v_))) .* (tem3 ./ tem2 + (tem1 ./ tem2).^2);
    
    hatz = z_tem(1: M) + 1j * z_tem(M + 1: 2 * M);
    hatv = max(v_tem(1: M) + v_tem(M + 1: 2 * M), 1e-10);
    %hatv = v_tem(1: M) + v_tem(M + 1: 2 * M);

    %v_tem = max(v_tem, 1e-10);
    hat_s = (hatz - Z) ./ (V);
    hat_tau = (V - hatv) ./ (V.^2);
end
function [umean, uvar] = estim(xo, v, wvar)
            % Get prior
            %umean0 = obj.mean0;
	        %uvar0 = max(eps,obj.var0); % avoid zero variances! 
            %M = obj.M ;
            
            % Compute posterior mean and variance
            AllSymbol = xo ;
            [size_n,size_m] = size(v) ; 
            v_temp = reshape(v, size_m*size_n, 1) ;
            wvar_temp = reshape(wvar, size_m*size_n, 1) ;
                        
             logpxr = bsxfun(@times, -1./(2*wvar), abs(bsxfun(@minus, v, AllSymbol)).^2) ;
             logpxr = bsxfun(@minus, logpxr, max(logpxr, [], 2) );
             pxr = exp(logpxr);
             pxr = bsxfun(@rdivide, pxr, sum(pxr,2) ) ; 
             umean = sum(bsxfun(@times, pxr, AllSymbol), 2) ; 
             uvar = sum(pxr .* abs(bsxfun(@minus, umean, AllSymbol)).^2, 2) ;
            
            % The above expressions are easily to understand 
            % speed up the peformance ---- 
            %if (~obj.maxSumVal)
                %logpxr = bsxfun(@times, -1./(2*wvar_temp), abs(bsxfun(@minus, v_temp, AllSymbol)).^2) ;
                % zhc: ????? logpxr = bsxfun(@times, -1./(wvar_temp), abs(bsxfun(@minus, v_temp, AllSymbol)).^2) ;
                %pxr = exp(logpxr);
                %pxrSum = sum(pxr,2) ;
               % umean = sum(bsxfun(@times, pxr, AllSymbol), 2)./pxrSum ;
                %uvar = sum(pxr .* abs(bsxfun(@minus, umean, AllSymbol)).^2, 2)./pxrSum ;
            %else
                %logpxr = bsxfun(@times, -1./(2*wvar_temp), abs(bsxfun(@minus, v_temp, AllSymbol)).^2) ;
                % zhc: ????  logpxr = bsxfun(@times, -1./(wvar_temp), abs(bsxfun(@minus, v_temp, AllSymbol)).^2) ;
                %[xxx, demoIdx] = max(logpxr, [], 2) ;
               % umean_temp = AllSymbol(demoIdx) ;
               % uvar_temp = wvar;
            %end
           
            %umean = reshape(umean_temp, size_n, size_m) ;
           % uvar = reshape(uvar_temp, size_n, size_m) ;  
            % val = reshape(val, size_n, size_m) ;  
        end