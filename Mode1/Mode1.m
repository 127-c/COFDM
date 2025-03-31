clc
clear all
%定义DAB mode 1 参数    %
%   帧持续时间：96ms   每帧包含4个CIF(不包含FIC)；每CIF有55296 Bit    %
%   FIC总数据率：96/kbps FIB：12个  MSC总数据率：2304/kbps    %
%   符号有效期： 1ms 副载波数目：1536   %
%   脚本只包含DQPSK 和 OFDM 调制 生成的DAB数据流链路，未进行卷积编码

% --- DAB模式一参数 ---
Con_off = 1;  % 不使用卷积编码
Con_on  = 0;
% 保护间隔点数（Δ=504T=246μs）
% 基础时间单位T=1/2048000秒
% 采样率=2.048 MHz
required_bits = 233472;%补上数据，不使用卷积编码不匹配
%参数设置
if Con_off  %不使用卷积编码
frame_duration = 96e-3; %帧持续时间
bitsrate_MSC = 2304e3;%MSC数据率
bitsrate_FIC = 96e3;%FIC数据率
T =1/2048000; 
fs = 1/T;


%生成 原始 比特流  
FIC_bits = randi([0,1], 96e3 * frame_duration, 1);      % 固定FIC=9,216 bits 

MSC_required = required_bits - length(FIC_bits);        % MSC需补足至233,472-9,216=224,256 bits 

MSC_bits = randi([0,1], MSC_required, 1);               % 修正MSC比特数 

Total_bits = [FIC_bits; MSC_bits];  %233,472 bits

end

if Con_on %使用卷积编码 
end

%------------ DQPSK 调制 -------------------

num_Subc = 1536;%有效子载波数

num_Loop = 76;%一共76个OFDM符号+1个0符号

SymQPSKtmp = reshape(Total_bits,2,num_Subc*num_Loop)'; 

%转为0-3；用于QPSK符号索引

SymQPSKtmptmp = bi2de(SymQPSKtmp,2,'left-msb');

%   --- 定义π/4 DQPSK查找表 ---
%00->pi/4 01->3*pi/4  10->-3*pi/4  11->-pi/4 
delta_phase_table = [pi/4 , 3*pi/4 , -3*pi/4 , -pi/4];

initial_ref_phase = 0;%初始化基准相位
current_ref_phase = initial_ref_phase;%前相位定义



%---生成 差分编码复数符号
diff_symbols = zeros( num_Loop*num_Subc , 1);%储存所有QPSK符号

for i=1 :num_Subc*num_Loop
    %获取输入索引
    idx = SymQPSKtmptmp(i);

    %获取相位变化
    delat_phase = delta_phase_table(idx + 1);

    %当前相位 = 变化相位 + 前相位
    current_phase = delat_phase + current_ref_phase;
switch SymQPSKtmptmp(i)
    case 0 
        amplitude = 0.8;
    case 1 
        amplitude = 1.0;
    case 2 
        amplitude = 1.2;
    case 3 
        amplitude = 1.4;
end 
    % 将相位转为复数（如π/4 → (1+1i)/√2）
    diff_symbols(i) =  amplitude *exp(1i *current_phase);
    
    %状态传递
    current_ref_phase = current_phase;

end

%-----OFDM调制-----------
%现在将DQPSK符号分配到有效子载波上
ifft_size = 2048;
num_subc_total = 1536;
cp_length = 504; %保护间隔长度

% 初始化频域符号矩阵（补零后的IFFT输入）
ofdm_symbols_freq = zeros(ifft_size , num_Loop); %每行代表一个OFDM符号频域数据

%计算1536载波在2048里的位置，然后前后都补上0
start_idx_ifft = (ifft_size - num_subc_total)/2 + 1 ;
end_idx_ifft = start_idx_ifft + num_subc_total-1;

for i =1 :num_Loop
    %循环处理每个OFDM符号
    start_data = (i-1)*num_Subc + 1;%第i个OFDM符号的1载波
    end_data = num_Subc*i;%第i个OFDM符号的第1536个载波
    amplitude_scaling = 1; % 示例缩放因子 
    data_symbols = diff_symbols(start_data : end_data)* amplitude_scaling;;%存储第i个OFDM符号的数据

    % 将1536数据补零到2048点（居中）
    ofdm_symbols_freq(start_idx_ifft:end_idx_ifft,  i) = data_symbols;
end

%-------IFFT与保护间隔（GI）添加
tx_signal = [];
for i = 1:num_Loop
    %对每个OFDM进行IFFT
    %频域 → 时域：将多个子载波的频域数据合并成一个时域波形
    ifft_data = ifft(ofdm_symbols_freq(:,i) , ifft_size);
    %保留2048个采样点
    ifft_data_valid = ifft_data(1:ifft_size);

    %CP
    cp = ifft_data_valid(end - cp_length +1 : end);%取末尾504个点作为保护间隔
    % 生成OFDM符号时域信号（504+2048=2552点）
    ofdm_symbol_time = [cp; ifft_data_valid];
    tx_signal = [tx_signal; ofdm_symbol_time];
end

% --- 生成零符号的时域信号（全零）---
null_symbol_time = zeros(2656, 1);  % 零符号持续时间为2656T

%0符号和OFDM符号拼接
tx_frame = [null_symbol_time; tx_signal];  % 总点数=2656+76×2552=196,608点

%-----绘制时域图像
t_total = (0:length(tx_frame) -1)/fs;%计算时间
% 绘制实部和虚部分量
figure;
subplot(2,1,1);%图像两行1列的第一个
plot(t_total, real(tx_frame));
xlabel('时间 (秒)');
ylabel('实部');
title('OFDM时域信号（实部）');
grid on;%网格线

subplot(2,1,2);
plot(t_total, imag(tx_frame));
xlabel('时间 (秒)');
ylabel('虚部');
title('OFDM时域信号（虚部）');
grid on;


% 提取第10个OFDM符号的数据 
symbol_index = 10; % 设置为10，提取第10个OFDM符号 
symbols_length = cp_length + ifft_size; % 每个OFDM符号的长度（CP + IFFT长度）
 
% 计算第10个符号的起始和结束索引 
start_idx = (symbol_index - 1) * symbols_length + 1;
end_idx = symbol_index * symbols_length;
 
% 提取第10个符号的数据 
tx_symbol_10 = tx_signal(start_idx:end_idx);
 
% 计算时间向量（仅针对第10个符号）
t_symbol = (0:symbols_length - 1) / fs; % 时间轴 
 
% 绘制实部和虚部分量 
figure;
subplot(2,1,1); % 上半部分绘制实部 
plot(t_symbol, real(tx_symbol_10));
xlabel('时间 (秒)');
ylabel('实部');
title('第10个OFDM符号时域信号（实部）');
grid on;
 
subplot(2,1,2); % 下半部分绘制虚部 
plot(t_symbol, imag(tx_symbol_10));
xlabel('时间 (秒)');
ylabel('虚部');
title('第10个OFDM符号时域信号（虚部）');
grid on;


%------------------接收端代码
%直接使用基带信号进行IQ调制，然后对有效符号做FFT
% --- 接收端处理（直接使用基带信号，无需下变频）---
% 假设接收到的信号是tx_frame（理想无噪声）

%----1.同步：检测0符号的起始位置
window_size  = 2656; %0符号长度

energy = movmean(abs(tx_frame).^2,window_size);%滑动能量检测
[~,sync_pos] = min(energy);%找到能量最低点

%-----2.截取有效信号（跳过0符号）
frame_start = sync_pos + window_size;%有效数据起始
%这里提取了除0符号以外的所有符号（相位同步符号也被提取进去了）%
rx_data = tx_frame(frame_start:end);%提取OFDM符号有效部分

%-----3.分割符号移除循环前缀
symbols_length = cp_length + ifft_size;
rx_symbols = reshape(rx_data , symbols_length , num_Loop).';%分割成OFDM符号

rx_symbols_valid = rx_symbols(:,cp_length+1:end);%去除CP

%-----4.FFT
rx_freq_symbols = zeros(num_Loop,ifft_size);%存储FFT结果

for i=1 :num_Loop
     %rx_freq_symbols( i, : )存储在i行所有列中
    rx_freq_symbols( i, : ) = fft(rx_symbols_valid( i, : ),ifft_size);%FFT-2048
end

%-----5.提取有效子载波（1536）
start_idx = (ifft_size - num_Subc)/2 + 1; % 有效子载波起始位置257

end_idx = start_idx + num_Subc - 1;       % 结束位置1792

rx_data_symbols = rx_freq_symbols(: , start_idx:end_idx); % 提取有效数据

symbol_idx = 10;
% 接收端FFT后的频域数据 
rx_freq_after_fft = rx_freq_symbols(symbol_idx, :).';
 
% 解调时考虑幅度变化 
rx_amplitude = abs(rx_freq_after_fft);
rx_phase = angle(rx_freq_after_fft);
 
% 根据幅度恢复原始数据 
% 示例：假设发送端使用了不同的幅度值 、

recovered_data = zeros(size(rx_amplitude));
for i = 1:length(rx_amplitude)
    if rx_amplitude(i) > 1.2 
        recovered_data(i) = 3; % 示例：最高幅度对应3 
    elseif rx_amplitude(i) > 1.0 
        recovered_data(i) = 2; % 示例：次高幅度对应2 
    elseif rx_amplitude(i) > 0.8 
        recovered_data(i) = 1; % 示例：中等幅度对应1 
    else 
        recovered_data(i) = 0; % 示例：最低幅度对应0 
    end 
end 
%6.绘制图像






% 在发射端添加IFFT之前的频谱图绘制 
% 选择一个OFDM符号进行分析 
symbol_index = 10;
freq_domain_tx = ofdm_symbols_freq(:, symbol_index);
% 计算幅度并进行归一化 
magnitude_tx = abs(freq_domain_tx) / max(abs(freq_domain_tx));
 
% 生成频率轴 
n_fft = ifft_size;
df = fs / n_fft;
f = (-n_fft/2:n_fft/2 - 1) * df;
 
% 绘制发射端频谱图 
figure;
plot(f, magnitude_tx);
xlabel('频率 (Hz)');
ylabel('归一化幅度');
title('IFFT之前的频谱图');
grid on;
 
% 在接收端添加FFT之后的频谱图绘制 
% 选择一个OFDM符号进行分析 
rx_freq_symbol = rx_freq_symbols(symbol_index, :);
% 计算幅度并进行归一化 
magnitude_rx = abs(rx_freq_symbol) / max(abs(rx_freq_symbol));
 
% 绘制接收端频谱图 
figure;
plot(f, magnitude_rx);
xlabel('频率 (Hz)');
ylabel('归一化幅度');
title('FFT之后的频谱图');
grid on;
