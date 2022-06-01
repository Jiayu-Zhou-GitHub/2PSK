function [output] = decision(o_NRZ)
%DECISION 此处显示有关此函数的摘要
%   此处显示详细说明
output=zeros(1,length(o_NRZ));        %声明变量的空间
%循环抽样判决
for n=1:length(o_NRZ)
    if o_NRZ(n)>0
        output(n)=1;
    elseif o_NRZ(n)<0
        output(n)=0;
    end
end
end

