function [MSE_error] = GAMP_Iter(obj, Input)
    eps = 1e-30;
    iter_num = Input.IterNum;
    %mes = Input.mes;

    N = Input.N;
    M = Input.M;
    H = obj.H;
    xo = obj.xo;

    %% array initialization
    hat_x = zeros(N, 1);
    hat_v = ones(N, 1);
    hat_s = zeros(M, 1);
    sqrH = abs(H).^2;
    sqrHt = sqrH';
    Ht = H';
    %V_old = zeros(M, 1);
    %Z_old = zeros(M, 1);


    MSE_error = zeros(iter_num, 1);

    hat_x_old = 0;
    hat_v_old = 0;

    for ii = 1: iter_num
        %% Output Nodes
        V = sqrH * hat_v;
        Z = H * hat_x - hat_s .* V;

        %[Z, Z_old] = damping(Z, Z_old, 0.8);
        %[V, V_old] = damping(V, V_old, 0.8);

        [hat_s, hat_tau] = estimator_2(Z, V, obj, Input);
        %% Input Nodes
        Sigma = (sqrHt * hat_tau).^(-1);
        R = hat_x + Sigma .* (Ht * hat_s);

        [hat_x, hat_v] = Estimator_1(xo, R, Sigma);
        if sum(isnan(hat_v)) > 0 
            hat_x = hat_x_old;
            hat_v = hat_v_old;
            break;
        end
        %hat_v = max(hat_v, eps);
        [hat_x, hat_x_old] = damping(hat_x, hat_x_old, 0.8);
        [hat_v, hat_v_old] = damping(hat_v, hat_v_old, 0.8);
        %hat_x_old = hat_x;
        %hat_v_old = hat_v;
        MSE_error(ii, 1) = mean(abs(hat_x - obj.x).^2);
    end
    %save('hat_x.mat', 'hat_x');

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
    log_posterior = bsxfun(@times, -1 ./ check_v, abs(bsxfun(@minus, xo, check_m).^2));
    log_posterior = bsxfun(@minus, log_posterior, max(log_posterior));  %防止溢出

    posterior = exp(log_posterior); 
    posterior = bsxfun(@rdivide, posterior, (sum(posterior, 2) + eps));   %得到标准PDF
    m = sum(bsxfun(@times, posterior, xo), 2);                    %计算PDF的均�?
    v = sum(posterior .* abs(bsxfun(@minus, m, xo).^2), 2);       %计算PDF的方�?
end

function [hat_s, hat_tau] = estimator_2(Z, V, obj, Input)
    %eps = 1e-30;
    nuw = Input.nuw;
    y = obj.y2;
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