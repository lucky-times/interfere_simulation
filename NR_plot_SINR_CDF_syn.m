function NR_plot_SINR_CDF_syn(SYS_config,UEs,activate_UE_id)
% ����SINR cdf����
% ��������� 
% SYS_config������
% UEs��UE�б�
% activate_UE_id����Ҫͳ�Ƶ�UE��id

%һЩ��ʼ��
all_SINR = [UEs.ul_dl_wideband_SINR];
ul_all_SINR = [UEs.dl_ul_wideband_SINR];

%% ���и�������
UE_SINR = all_SINR(activate_UE_id);
UE_SINR = UE_SINR(UE_SINR>-inf);
[range,output] = NR_plot_CDF(UE_SINR);
if SYS_config.photo
    figure;
    subplot(1,2,1);
    plot(range,output,'b');
    grid on;
    xlabel('SINR��dB��');
    ylabel('F(x)');
    title('UL_to_DL SINR CDF');
end
%% �������и�������SINR����
UE_SINR = sort(UE_SINR);
index = 1:length(UE_SINR)/100:length(UE_SINR);
index = round(index);
R_sinr = UE_SINR(index);
R_sinr = [R_sinr max(UE_SINR)];
fid=fopen('.\calibration\R_sinr_syn.txt','wt');
fprintf(fid,'%f\n',R_sinr);
fclose(fid);

%% ���и�������
ul_UE_SINR = ul_all_SINR(activate_UE_id);
ul_UE_SINR = ul_UE_SINR(ul_UE_SINR>-inf);
[range,output] = NR_plot_CDF(ul_UE_SINR);
if SYS_config.photo
    subplot(1,2,2);
    plot(range,output);
    xlabel('SINR��dB��');
    ylabel('F(x)');
    title('DL_to_UL SINR CDF');
end

ul_UE_SINR = sort(ul_UE_SINR);
index = 1:length(ul_UE_SINR)/100:length(ul_UE_SINR);
index = round(index);
R_sinr = ul_UE_SINR(index);
R_sinr = [R_sinr max(ul_UE_SINR)];
fid=fopen('.\calibration\R_sinr_ul_syn.txt','wt');
fprintf(fid,'%f\n',R_sinr);
fclose(fid);
end