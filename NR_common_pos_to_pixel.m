function [ pos_pixel,pos_pixel_exact] = NR_common_pos_to_pixel( pos,       roi_min, data_res)
% position��������ֵλ��ת��������λ�ã�����һ����Ӧ
% �����ص�λ��ת����λ�þ���ֵ
% ���������
% pos��ʵ�ʵ�λ��
% roi_min��ROI�����ĺ����²�ֵ
% data_res�����ݷֱ���
%
% ���������
% pos_pixel�����ص��λ��


pos_pixel(:,1) = floor((pos(:,1)-roi_min(1))/data_res)+1;
pos_pixel(:,2) = floor((pos(:,2)-roi_min(2))/data_res)+1;

pos_pixel_exact(:,1) = ((pos(:,1)-roi_min(1))/data_res)+1;
pos_pixel_exact(:,2) = ((pos(:,2)-roi_min(2))/data_res)+1;

