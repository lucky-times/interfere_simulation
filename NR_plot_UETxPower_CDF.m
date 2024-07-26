function NR_plot_UETxPower_CDF(SYS_config,UEs,activate_UE_id)
% ����UE���з��书�ʵ�cdf
% ���������
% SYS_config��LTE����ϵͳ
% UEs��UE�б�

for u_=1:length(activate_UE_id)
    UE_tx_power(:,u_) = UEs(activate_UE_id(u_)).ul_ad_TX_power; % ���빦�غ�ķ��书��
end
UE_tx_power = reshape(UE_tx_power,1,[]);
UE_tx_power = 10*log10(UE_tx_power*1000); % ����ֵתΪdBm

if SYS_config.photo
    [range,output] = NR_plot_CDF(UE_tx_power);
    figure;
    plot(range,output);
    xlabel('UL Tx power(dBm)');
    ylabel('F(x)');
    title('UL ad frequency Tx power CDF');
end
%% ���������ݱ���Ϊtxt�ĵ�
UE_tx_power = sort(UE_tx_power);
index = 1:length(UE_tx_power)/100:length(UE_tx_power);
index = round(index);
R_ue_power = UE_tx_power(index);
R_ue_power = [R_ue_power max(UE_tx_power)];
fid=fopen('.\calibration\R_ue_power.txt','wt');
fprintf(fid,'%f\n',R_ue_power);
fprintf('���ĵ���¼ͳ�Ƹ���UE���غ�ķ��书��\n');
fclose(fid);

end