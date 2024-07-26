function NR_plot_IoT_CDF(SYS_config,IoT)
% ����IoTcdf����
% ���������
% SYS_config������ϵͳ
% IoT�����Ŷ������ı�ֵ

if SYS_config.photo
    % ����IoT cdf����
    % ��ϵͳ��ToT�еĸ��Ž�ΪͬƵ����
    % ˫ϵͳ�еĸ��Ű���ͬƵ����Ƶ����
    IoT_1 = IoT;
    [range,output] = NR_plot_CDF(IoT_1);
    figure;
    plot(range,output,'b');
    grid on;
    xlabel('IoT��dB��');
    ylabel('F(x)');
    title('IoT CDF');
    % ��������Ϊ�˽���ͬ�������жȵĽ������ͬһ��fig��
    file_name1 = ['.\tmp\','IOT_range_',num2str((SYS_config.PCe_param)^2)];% ������tmp�У������Ҫ����ͬ���ش�������µ�ͼ�Σ�����д���ű�����
    file_name2 = ['.\tmp\','IOT_output_',num2str((SYS_config.PCe_param)^2)];
    save(file_name1,'range');
    save(file_name2,'output');
end

end