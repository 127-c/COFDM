clc
clear all
%定义DAB mode 1 参数    %
%   帧持续时间：96ms   每帧包含4个CIF(不包含FIC)；每CIF有55296 Bit    %
%   FIC总数据率：96/kbps FIB：12个  MSC总数据率：2304/kbps    %
%   符号有效期： 1ms 副载波数目：1536   %
%   脚本只包含DQPSK 和 OFDM 调制 生成的DAB数据流链路，未进行卷积编码以及符号同步
%   总子载波1536中，仅1200个用于数据，其余用于同步和保护间隔。
Con_off = 1;
Con_on  = 0;

%参数设置
if Con_off  %不使用卷积编码
frame_duration = 96e-3; %帧持续时间
num_CIF_frame = 4;%CIF个数
num_FIB_frame = 12;%FIB个数
bitsrate_FIB = 96e3;%FIB数据率
bits_per_CIF = 55296;%每CIF 包含Bits数
bits_per_FIB = (bitsrate_FIB*frame_duration)/12; %每FIB 包含Bits数
bits_per_symbol = 2400;%每符号bit数 2304 bits+96 bits 1200*2bits/子载波数



%生成 原始 比特流  
MSC_bits = randi([0,1] , num_CIF_frame*bits_per_CIF ,1 );%MSC

FIC_bits = randi([0,1], bits_per_FIB*num_FIB_frame , 1 );%FIC

end

if Con_on %使用卷积编码 
end

%------------ DQPSK 调制 -------------------

Total_bits = [FIC_bits;MSC_bits]';% 230400bits

num_Subc = 1200;%有效子载波数

num_Loop = length(Total_bits) / (num_Subc*2);%一共96个OFDM符号

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
valid_subc = 1:1200;%有效载波索引
cp_length = 384; %保护间隔长度

% 初始化频域符号矩阵（补零后的IFFT输入）
ofdm_symbols_freq = zeros(ifft_size , num_Loop); %每行代表一个OFDM符号频域数据

%计算1536载波在2048里的位置，然后前后都补上0
start_idx_ifft = (ifft_size - num_subc_total)/2 + 1 ;
end_idx_ifft = start_idx_ifft + num_subc_total-1;

%计算1200在1536载波里的位置
start_idx_total = (num_subc_total - num_Subc)/2 + 1;
end_idx_total = start_idx_total + num_Subc - 1;
for i =1 :num_Loop
    %循环处理每个OFDM符号
    start_data = (i-1)*num_Subc + 1;%第i个OFDM符号的1载波
    end_data = num_Subc*i;%第i个OFDM符号的第1200个载波

    data_symbols = diff_symbols(start_data : end_data);%存储第i个OFDM符号的数据

    % 在总子载波中居中填充有效数据（1536长度，中间1200为数据，其余补零）
    %填充第i个OFDM符号的中间1200个载波的数据
    total_subc_data = zeros(num_subc_total, 1 );%初始1536个数据

    %填充中间1200到1536
    total_subc_data(start_idx_total : end_idx_total) = data_symbols;

    % 将1536数据补零到2048点（居中）
    ofdm_symbols_freq(start_idx_ifft:end_idx_ifft,  i) = total_subc_data;
end

%-------IFFT与保护间隔（GI）添加
tx_signal = [];
for i = 1:num_Loop
    %对每个OFDM进行IFFT
    %频域 → 时域：将多个子载波的频域数据合并成一个时域波形
    ifft_data = ifft(ofdm_symbols_freq(:,i) , ifft_size);
    %取前1536个采样点
    ifft_data_valid = ifft_data(1:num_subc_total);

    %CP
    cp = ifft_data_valid(end - cp_length : end);%取末尾384个点作为保护间隔
    
    tx_signal = [tx_signal ; cp ; ifft_data_valid];
end
