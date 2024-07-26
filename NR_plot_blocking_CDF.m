function NR_plot_blocking_CDF( SYS_config,interference_level )
% ���ƴ������� cdf����
% ���������
% SYS_config��LTE����ϵͳ
% interference_level�� �����ź�

if SYS_config.photo
    if SYS_config.isDouble
        % ���������ź�CDF
        interference_level_1 = interference_level;
        interference_level_1 = 10.*log10(interference_level_1.*1000);
        [range,output] = NR_plot_CDF(interference_level_1);
        
        % Ϊ����ȡ��99.99%�ĵ㣬�����߽���ƽ��
        ft = fittype( 'smoothingspline' );
        opts = fitoptions( ft );
        opts.SmoothingParam = 1;
        [fitresult, ~] = fit( range', output', ft, opts );
        range = linspace(min(range),max(range),length(range)*100);
        output = feval(fitresult,range);
        
        figure;
        plot(range,output,'b');
        axis([-inf,inf,0,1]);
        grid on;
        xlabel('Blocking��dBm��');
        ylabel('F(x)');
        title('in-band blocking signal CDF');
        % ��������Ϊ�˽���ͬ�������жȵĽ������ͬһ��fig��
        file_name1 = ['.\tmp\','Blocking_range_',num2str((SYS_config.PCe_param)^2)];% ������tmp�У������Ҫ����ͬ���ش�������µ�ͼ�Σ�����д���ű�����
        file_name2 = ['.\tmp\','Blocking_output_',num2str((SYS_config.PCe_param)^2)];
        save(file_name1,'range');
        save(file_name2,'output');
    end
end


end

