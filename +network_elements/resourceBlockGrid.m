classdef resourceBlockGrid < handle
    % RB�࣬��Ҫ�����Ǹ�ÿ��RB���书��
    % ����ֻ���õ�������Ϊһ�ֶ���ʹ�ô���ķ�����
    % ���ǿ��ǵ��������Ǹ�ÿ��UEȫ����RB����������ACIRģ����LTE��һ���Ļ���
    % ������һ�������б�Ҫ��
    
    properties
        % BS��UE�յĹ��ʷ��䣬��λw
        power_allocation
        % UE��BS�յĹ��ʷ��䣬��λw
        ue_power_allocation
        % ��Ƶ��
        power_allocation_signaling
        ue_power_allocation_signaling
        % RB��
        n_RB
    end
    
    methods
        
        % ���캯��
        function obj = resourceBlockGrid(n_RB)
            obj.power_allocation = zeros(n_RB,1); 
            obj.power_allocation_signaling = zeros(n_RB,1);
            obj.ue_power_allocation = zeros(n_RB,1);
            obj.ue_power_allocation_signaling = zeros(n_RB,1);
            obj.n_RB = n_RB;
        end
        
        % ÿ��RB���䵽�Ĺ�����һ����
        function set_homogeneous_power_allocation(obj,SYS_config,power_in_watts_data,power_in_watts_signaling)
            % power_in_watts_data = eNodeBs(b_).sectors(s_).max_power
            % power_in_watts_signaling = eNodeBs(b_).sectors(s_).signaling_power
            
            % ��վ��
            power_per_RB_data                 = power_in_watts_data / obj.n_RB;%����ƽ�������ÿ��RB
            power_per_RB_signaling            = power_in_watts_signaling / obj.n_RB;
            obj.power_allocation(:)           = power_per_RB_data;
            obj.power_allocation_signaling(:) = power_per_RB_signaling;
            
            UE_tx_power = SYS_config.UE_tx_power;
            ue_max_power = UE_tx_power * (1-SYS_config.signaling_ratio);
            ue_signaling_power = UE_tx_power * SYS_config.signaling_ratio;
            
            % UE��
            % ʵ���ϻ�վ��UEֻ�ǹ��ʲ�ͬ
            ue_power_per_RB_data = ue_max_power / obj.n_RB;%����ƽ�������ÿ��RB
            ue_power_per_RB_signaling = ue_signaling_power / obj.n_RB;
            obj.ue_power_allocation(:)           = ue_power_per_RB_data;
            obj.ue_power_allocation_signaling(:) = ue_power_per_RB_signaling;
        end

        % Prints some info about this object
        function print(obj)
            fprintf('n_RB=%d\n',obj.n_RB);
        end
        
    end
end
