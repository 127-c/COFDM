clc
clear all
%定义DAB mode 1 参数    %
%   帧持续时间：96ms   每帧包含4个CIF(不包含FIC)；每CIF有55296 Bit    %
%   FIC总数据率：96/kbps FIB：12个  MSC总数据率：2304/kbps    %
%   符号有效期： 1ms 副载波数目：1536   %
%   脚本只包含DQPSK 和 OFDM 调制 生成的DAB数据流链路，未进行卷积编码以及符号同步
%   总子载波1536中，仅1200个用于数据，其余用于同步和保护间隔。

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

    % 将相位转为复数（如π/4 → (1+1i)/√2）
    diff_symbols(i) = exp(1i *current_phase);
    
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

    data_symbols = diff_symbols(start_data : end_data);%存储第i个OFDM符号的数据

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
tx_frame = [null_symbol_time; tx_signal];  % 总点数=2656 + 76×2552=196,608点


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
