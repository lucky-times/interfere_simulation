function [ pos ]                      = NR_common_pixel_to_pos( pos_pixel, roi_min, data_res)
% �����ص�λ��ת����λ�þ���ֵ
% ���������
% pos_pixel�����ص��λ��
% roi_min��ROI�����ĺ����²�ֵ
% data_res�����ݷֱ���
%
% ���������
% pos��ʵ�ʵ�λ��

pos(:,1) = (pos_pixel(:,1)-1)*data_res+roi_min(1);
pos(:,2) = (pos_pixel(:,2)-1)*data_res+roi_min(2);
