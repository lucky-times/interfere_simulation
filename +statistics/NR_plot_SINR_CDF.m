function NR_plot_SINR_CDF(SYS_config,UEs,activate_UE_id)
% ����SINR cdf����
% ��������� 
% SYS_config������
% UEs��UE�б�
% activate_UE_id����Ҫͳ�Ƶ�UE��id

%һЩ��ʼ��
all_SINR = [UEs.wideband_SINR];
ul_all_SINR = [UEs.ul_wideband_SINR];

%% ����
UE_SINR = all_SINR(activate_UE_id);
[range,output] = NR_plot_CDF(UE_SINR);
if SYS_config.photo
    figure;
    subplot(1,2,1);
    plot(range,output,'b');
    grid on;
    xlabel('SINR��dB��');
    ylabel('F(x)');
    title('DL SINR CDF');
end
%% ��������SINR����
UE_SINR = sort(UE_SINR);
index = 1:length(UE_SINR)/100:length(UE_SINR);
index = round(index);
R_sinr = UE_SINR(index);
R_sinr = [R_sinr max(UE_SINR)];
fid=fopen('.\calibration\R_sinr.txt','wt');
fprintf(fid,'%f\n',R_sinr);
fclose(fid);

%% ����
ul_UE_SINR = ul_all_SINR(activate_UE_id);
[range,output] = NR_plot_CDF(ul_UE_SINR);
if SYS_config.photo
    subplot(1,2,2);
    plot(range,output);
    xlabel('SINR��dB��');
    ylabel('F(x)');
    title('UL SINR CDF');
end

ul_UE_SINR = sort(ul_UE_SINR);
index = 1:length(ul_UE_SINR)/100:length(ul_UE_SINR);
index = round(index);
R_sinr = ul_UE_SINR(index);
R_sinr = [R_sinr max(ul_UE_SINR)];
fid=fopen('.\calibration\R_sinr_ul.txt','wt');
fprintf(fid,'%f\n',R_sinr);
fclose(fid);
end