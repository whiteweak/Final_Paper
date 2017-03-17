% 2.2.2 将二进制编码转化为十进制数(2)
% decodechrom.m函数的功能是将染色体(或二进制编码)转换为十进制，参数spoint表示待解码的二进制串的起始位置
% (对于多个变量而言，如有两个变量，采用20为表示，每个变量10为，则第一个变量从1开始，另一个变量从11开始。本例为1)，
% 参数1ength表示所截取的长度（本例为10）。
%遗传算法子程序
%Name: decodechrom.m
%将二进制编码转换成十进制

function pop2=decodechrom(pop,spoint,length)          %这里的spoint的值为 1，length的值为 10，这两个值的含义是:
                                                      %将pop数列中的每一行的第1个到第10个二进制数转换成相应的十进制数
                                                      %结果是20*1的矩阵
pop1=pop(:,spoint:spoint+length-1);                   %pop1等价于pop，这里的目的是为了不改变原始pop的值
pop2=decodebinary(pop1);                              %pop2表示的则是一个20x1的矩阵，为pop的元素转换成10进制之后的表现形式