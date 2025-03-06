
function [Trellis, Coded_bits, Bits] = COFDM_Code(N, Rate, Block_Num, State)  
    Bits = randi([0 1], 1, N * Block_Num);  %生成原始的bit流
    Trellis = [];  % 初始化Trellis为空  
    
    if Rate == 1/2  
        % 约束长度3，生成多项式7(111),5(101)  
        Trellis = poly2trellis(3, [7, 5]);  %把多项式转化为网格表达式
    elseif Rate == 2/3  
        % 多输入卷积码（约束长度2,2）  
        Trellis = poly2trellis([2,2], [2 1 2; 1 3 3]);%相当于定义寄存器的结构
    else  %多输入的卷积码，有两个分支，输出状态数为2的（1+1）次方，1为每个分支的寄存器个数
        error('不支持的编码速率: Rate必须为1/2或2/3');  % 抛出错误而非静默返回  
        return
    end  
    
    % 计算实际编码输出长度  
    code_rate = 1/Rate;%计算冗余放大系数  %rate若为2/3 1.5 表示 1 bit输出1.5bit数据
    Coded_bits = zeros(1, length(Bits) * code_rate);  %此时，生成的0的个数为实际的码长
    
    % 分段编码（确保输入长度可被State整除）  
    block_size = State;  %output状态个数
    for a = 1:block_size:length(Bits)  %分段，步长为block_size(4)
        end_idx = min(a + block_size -1, length(Bits));  %防止索引越界
        encoded_block = convenc(Bits(a:a+State-1), Trellis);%分段编码  
        start = (a-1)*code_rate +1;  %输出位置计算
        Coded_bits(start : start + length(encoded_block)-1) = encoded_block;%拼接起来编码流  
    end  
     % 创建可视化窗口  
    fig_handle = figure('Name', '卷积编码数据可视化', 'Position', [100 100 1200 800]);  

    % 1. 原始与编码比特流对比（前200比特）  
    subplot(3,1,1);  
    stem(Bits(1:200), 'filled', 'MarkerSize',4, 'LineWidth',1.5);  
    hold on;  
    stem(Coded_bits(1:200*code_rate), '^', 'LineWidth',1.2, 'Color',[0.8 0.2 0.2]);  
    title(['比特流对比 | 码率 ' num2str(Rate)]);  
    xlabel('比特序号'); ylabel('比特值');  
    legend('原始比特', '编码后比特', 'Location','northeast');  
    grid on;   
end  
