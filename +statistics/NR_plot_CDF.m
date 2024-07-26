function [range,output] = NR_plot_CDF(input)
% ��Ҫ��CDF�����ݽ���Ԥ����
% �������
% input������������
% �������
% range��CDFͼ�е�X��
% output��CDFͼ�е�Y��

min_ = floor(min(input));
max_ = ceil(max(input));
range = min_:1:max_;
h = hist(input,range);
output = cumsum(h)/sum(h);
end