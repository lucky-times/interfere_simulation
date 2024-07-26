function NR_plot_receive_power_CDF( SYS_config,UE_Rx_power,UE_Rx_intf_power,BS_Rx_power,BS_Rx_intf_power )
% ͳ�ƽ��յ�Ĺ���
% �������
% SYS_config����������
% UE_Rx_power��UE�˽��չ���
% UE_Rx_intf_power��UE�˽��յĸ��Ź���
% BS_Rx_power��BS�˽��չ���
% BS_Rx_intf_power��BS�˽��յĸ��Ź���

UE_Rx_power = 10.*log10(UE_Rx_power);
UE_Rx_intf_power = 10.*log10(UE_Rx_intf_power);
BS_Rx_power = 10.*log10(BS_Rx_power);
BS_Rx_intf_power = 10.*log10(BS_Rx_intf_power);
% plot
if SYS_config.photo
    [range,output] = NR_plot_CDF(UE_Rx_power);
    figure;
    subplot(2,2,1);
    plot(range,output,'b');
    grid on;
    xlabel('UE Rx power��dBw��');
    ylabel('F(x)');
    
    [range,output] = NR_plot_CDF(UE_Rx_intf_power);
    subplot(2,2,2);
    plot(range,output,'b');
    grid on;
    xlabel('UE Rx intf power��dBw��');
    ylabel('F(x)');
    
    try % InH�����£�UE�������غ�ʹ�û�վ���չ���һ����������CDF����
        [range,output] = NR_plot_CDF(BS_Rx_power);
        subplot(2,2,3);
        plot(range,output,'b');
        grid on;
        xlabel('BS Rx power��dBw��');
        ylabel('F(x)');
    catch
        
    end
    
    [range,output] = NR_plot_CDF(BS_Rx_intf_power);
    subplot(2,2,4);
    plot(range,output,'b');
    grid on;
    xlabel('BS Rx intf power��dBw��');
    ylabel('F(x)');
    
end

% svae
save_flie(UE_Rx_power,'R_UE_Rx_power');
save_flie(UE_Rx_intf_power,'R_UE_Rx_intf_power');
save_flie(BS_Rx_power,'R_BS_Rx_power');
save_flie(BS_Rx_intf_power,'R_BS_Rx_intf_power');

end

function save_flie(data,file_name)
data = sort(data);
index = 1:length(data)/100:length(data);
index = round(index);
R_data = data(index);
R_data = [R_data max(data)];
path = ['.\calibration\',file_name,'.txt'];
fid=fopen(path,'wt');
fprintf(fid,'%f\n',R_data);
fclose(fid);
end

