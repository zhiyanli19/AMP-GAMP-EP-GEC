function obj = EP_MIMO_System(Input)
	
	mod_size = Input.mod_size;
	N = Input.N;
	M = Input.M;
	%bit = Input.bit;
	nuw = Input.nuw;

	%% Generate xo
	obj.xo = Constellations(mod_size); % Generate original constellation sets

	%% Generate x
	sym = modem.qammod(2^mod_size);
    normal_scal = 1 / sqrt((2 / 3) * (2^mod_size - 1)) ;      %QAM normalization 
	sym.input = 'bit';                                  %the type of input data shoud be 'bit'
	sym.symbolorder = 'gray';
	information_bit = round(rand(N * mod_size, 1)) ;      %�������ж���������
	information_sym = modulate(sym, information_bit);   %�Զ��������н���QAM����
	x = normal_scal * reshape(information_sym, N, 1);   %�õ���һ����QAM�����ź�

    %% Channel
    H=(randn(M,N)+1j*randn(M,N))/sqrt(2*N);

%% Noise
    w=sqrt(nuw/2)*(randn(M,1)+1j*randn(M,1));
	%% Channel
	%H = (randn(M, N) + 1j * randn(M, N)) / sqrt(double(2 * N));
	%% Noise
	%w = sqrt(double(nuw / 2)) * (randn(M, 1) + 1j * randn(M, 1));   %������˹����
    %w = normrnd(0, sqrt(nuw / 2), [M, 1]) + 1j * normrnd(0, sqrt(nuw / 2), [M, 1]);
    %w = normrnd(0, nuw, [M, 1]) + 1j * normrnd(0, nuw, [M, 1]);
	%% Uncoded system
	z = H * x;
	z_w = z + w;
	%[y, quan_step] = AMP_Quantization(z_w, bit);
    y = z_w;
	%% load parameters
	obj.x = x;
	obj.H = H;
	obj.y = y;
	obj.information_bit = information_bit;
	%obj.quan_step = quan_step;
end