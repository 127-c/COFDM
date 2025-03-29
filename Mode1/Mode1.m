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
num_data_symbols = 96;%每帧符号数（OFDM符号）
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



%数据符号生成
zero_symbol = zeros(num_Subc,1);%生成 0 符号  




