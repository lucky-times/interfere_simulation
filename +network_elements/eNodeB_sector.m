classdef eNodeB_sector < handle
    % С��ʵ��

    properties
        id                           % С����eNodeB���������
        eNodeB_id                    % ��sector������sector�е�id�����Ǳ�sector����enb��id
        parent_eNodeB                % С��������eNodeB
        azimuth                      % С�������Ȱ����ת��
        antenna                      % С������
        attached_UEs = 0;            % С�����ӵ��û���
        max_power                    % С�����ķ��书��
        signaling_power              % Power dedicated to signaling. Counted as always in use��Ϊ0
        scheduler                    % С���ĵ�����
        RB_grid                      % С������Դ������
        nTX                          % С���ģ����⣩������
        in_interf_eNodeB_sectors             % ��һ��ϵͳ�ĸ�������
        ad_interf_eNodeB_sectors                % �ڶ���ϵͳ�ĸ�������
        
        always_on = true;            % �Ƿ�С�������������Ƿ����û����룩
        
        last_received_feedback       % ��������źŷ���
        attached_UEs_vector          % ����UE���б�
        
        feedback_trace               % ����UE���б�
        sector_trace                 % ����UE���б�
        
        zero_delay_feedback = false; % ����UE���б�

        macroscopic_pathloss_model   % С��ʹ�õ�·��ģ��
        

        unquantized_CQI_feedback = false;
        
        % extra params that may not be always used (mainly extra information from the network planning tool)
        transmitter     % This points to the pathloss file that will actually be used
        frequency_band  % С��ʹ�õ�Ƶ��
        antenna_name    % С����������
        antenna_type    % С��ʹ�õ���������
        tx_height       % С�������߸߶�
        electrical_downtilt % С�����ߵĵ��������
        mechanical_downtilt % С�����ߵĻ�е�����
    end

    methods
        
        function print(obj)
            fprintf(' Sector %d: ',obj.id);
            fprintf('%s %ddB %d\n',obj.antenna.antenna_type,obj.antenna.mean_antenna_gain,obj.azimuth);
            fprintf('  ');
            obj.scheduler.print;
            fprintf('  UEs: ');
            if ~isempty(obj.users)
                current_last_node = obj.users(1);
                while ~isempty(current_last_node.Next)
                    % Do something
                    fprintf('%d ',current_last_node.Data.id);
                    current_last_node = current_last_node.Next;
                end
                % Do something for the last node
                fprintf('%d ',current_last_node.Data.id);
            end
            fprintf('\n');
            fprintf('  '); obj.RB_grid.print;
        end
        
        function clear(obj)
            obj.parent_eNodeB              = [];
            obj.antenna                    = [];
            obj.attached_UEs_vector        = [];
            obj.scheduler                  = [];
            obj.RB_grid                    = [];
            obj.in_interf_eNodeB_sectors           = [];
            obj.ad_interf_eNodeB_sectors = [];
            obj.last_received_feedback     = [];
            obj.feedback_trace             = [];
            obj.sector_trace               = [];
            obj.macroscopic_pathloss_model = [];
            obj.transmitter                = [];
            obj.frequency_band             = [];
        end
        
        %% ��ָ���û����뵽С����
        function attachUser(obj,user)
            if isempty(obj.attached_UEs_vector)
 
                obj.attached_UEs_vector  = user;
                obj.attached_UEs         = obj.attached_UEs + 1;
            else

                current_UEs = [obj.attached_UEs_vector.id];
                if ~sum(current_UEs==user.id)
                    obj.attached_UEs_vector = [obj.attached_UEs_vector user];
                    obj.attached_UEs        = obj.attached_UEs + 1;
                end
            end
            
            % UE������attach����
            user.attached_sector_idx = obj.id;
            user.attached_site       = obj.parent_eNodeB;
            user.attached_eNodeB     = obj;
        end
        
        %% ��ָ���û���С���з���
        function deattachUser(obj,user)
            if ~isempty(obj.attached_UEs_vector)
                current_UEs  = [obj.attached_UEs_vector.id];
                UE_idx       = (current_UEs==user.id);
                UE_in_eNodeB = sum(UE_idx);
                
                if UE_in_eNodeB>0
                    obj.attached_UEs_vector = obj.attached_UEs_vector(~UE_idx);
                end
                obj.attached_UEs = obj.attached_UEs - UE_in_eNodeB;
            end
        end
        
         %% �ж�ĳ���û��Ƿ�����С��
        function is_attached = userIsAttached(obj,user)
            if ~isempty(obj.users)
                current_UEs  = [obj.attached_UEs_vector.id];
                UE_idx       = (current_UEs==user.id);
                UE_in_eNodeB = sum(UE_idx);
                is_attached  = logical(UE_in_eNodeB);
            else
                is_attached = false;
            end
        end       
        
        
        function clear_non_basic_info(obj)
            obj.antenna                    = [];
            obj.attached_UEs               = [];
            obj.attached_UEs_vector        = [];
            obj.scheduler                  = [];
            obj.RB_grid                    = [];
            obj.last_received_feedback     = [];
            obj.feedback_trace             = [];
            obj.sector_trace               = [];
            obj.macroscopic_pathloss_model = [];
        end
    end
end
