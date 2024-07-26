function NR_plot_signal_and_intf_CL_CDF(signal_CL,interfering_CL,another_interfering_CL,activate_UE_id)
% ����ÿ��ϵͳCL��cdf�Լ������źŵ�CL
% ���������
% signal_CL�������źŵ�������
% interfering_CL��ͬƵ�����źŵ�CL
% another_interfering_CL����Ƶ�����źŵ�CL
% UEs��UE�б�

n_UE = length(activate_UE_id);
% ȷ������������ź�CL��������Ϊ�˻��101�����ݱ���У׼��
s_CL = -signal_CL';
s_CL = sort(s_CL);
index = 1:length(s_CL)/100:length(s_CL);
index = round(index);
% ��������ź�
R_s_CL = s_CL(index);
R_s_CL = [R_s_CL max(s_CL)];
fid=fopen('.\calibration\R_s_CL.txt','wt');
fprintf(fid,'%f\n',R_s_CL);
fclose(fid);
% ȷ��������ͬƵ�����ź�CL������
CL = reshape(interfering_CL,1,n_UE*size(interfering_CL,2));
CL = -CL;
CL = sort(CL);
index = 1:length(CL)/100:length(CL);
index = round(index);
% ����ͬƵ��CL
R_CL = CL(index);
R_CL = [R_CL max(CL)];
fid=fopen('.\calibration\R_CL.txt','wt');
fprintf(fid,'%f\n',R_CL);
fclose(fid);
% ȷ����������Ƶ�����ź�CL������
another_CL = reshape(another_interfering_CL,1,n_UE*size(another_interfering_CL,2));
another_CL = -another_CL;
another_CL = sort(another_CL);
index = 1:length(another_CL)/100:length(another_CL);
index = round(index);
% ������Ƶ��CL
R_another_CL = another_CL(index);
R_another_CL = [R_another_CL max(another_CL)];
fid=fopen('.\calibration\R_another_CL.txt','wt');
fprintf(fid,'%f\n',R_another_CL);
fclose(fid);
end