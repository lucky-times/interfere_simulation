classdef naming
    % ��������ļ�����
    
    properties
        % ʱ��
        the_date
        % ʱ��ת�����ַ���
        date_time_string
    end
    
    methods
        % ���캯��������ǰʱ��ת�����ַ�����������
        function obj = naming
            obj.the_date = clock;
            obj.date_time_string = sprintf('%04d%02d%02d_%02d%02d%02d',...
                obj.the_date(1),...                     % Date: year
                obj.the_date(2),...                     % Date: month
                obj.the_date(3),...                     % Date: day
                obj.the_date(4),...                     % Date: hour
                obj.the_date(5),...                     % Date: minutes
                floor(obj.the_date(6)));                % Date: seconds
        end
        
        % ���н���ļ�������
        function results_filename = results_file(obj,SYS_config)
            if strcmp(SYS_config.results_file,'auto')
                if SYS_config.frequency/1e9 >= 1
                    this_freq = sprintf('%3.2fGHz',SYS_config.frequency/1e9);
                else
                    this_freq = sprintf('%3.0fMHz',SYS_config.frequency/1e6);
                end
                
                this_bw = sprintf('%gfMHz',SYS_config.bandwidth/1e6);
                
                if SYS_config.keep_UEs_still
                    UE_speed = 0;
                else
                    UE_speed = SYS_config.UE_speed;
                end
                
                % Output file name
                results_filename = sprintf('%s_freq_%s_bw_%.1fKmph_%dTTIs_%s.mat',...
                    this_freq,...                             % Ƶ��
                    this_bw,...                               % ����
                    UE_speed*3.6,...                          % �ٶ�
                    SYS_config.simulation_time_tti,...        % TTI
                    obj.date_time_string);                    % ʱ��
            else
                results_filename = SYS_config.results_file;
            end
            
            % Check for .mat extension
            if ~strcmp(results_filename((end-3):(end)),'.mat')
                results_filename = [results_filename '.mat'];
            end
        end
        
        function beam_cache_filename = beam_cache_file(obj,SYS_config)
            
            if SYS_config.frequency/1e9 >= 1
                this_freq = sprintf('%3.2fGHz',SYS_config.frequency/1e9);
            else
                this_freq = sprintf('%3.0fMHz',SYS_config.frequency/1e6);
            end
            
            if isempty(SYS_config.ISD)
                ISD = 'null';
            else
                ISD = sprintf('%d',SYS_config.ISD);
            end
            
            beam_cache_filename = sprintf('%s_%s_ISD%s_UE%d_isDouble%d_ShiftMode%d_isHandover%d_AntennaMode%d_AttatchMode%d_isWraparound%d_isNewUMA%d_isManhattan%d_isFemto%d',...
               this_freq,...                    % Ƶ��
               SYS_config.scene_type,...        % ����
               ISD,...                          % ��վ�����
               SYS_config.UE_per_eNodeB,...     % ÿ��С��UE��
               SYS_config.isDouble,...          % �Ƿ�˫ϵͳ
               SYS_config.shift_mode,...        % ƫ������
               SYS_config.hand_over,...         % hand over
               SYS_config.antenna_mode,...      % ����ģʽ
               SYS_config.attatch_mode,...      % ����ģʽ
               SYS_config.isWraparound,...      % wraparound
               SYS_config.isNewUMA,...          % ������UE����ѡ��
               SYS_config.isManhattan,...       % ������
               SYS_config.isFemto);             % Femto
        end
        
    end
end

