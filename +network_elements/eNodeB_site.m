classdef eNodeB_site < handle
    % eNodeBʵ�壨��վ��
    
    properties
        id                  % eNodeB ��id
        pos                 % eNodeB��λ�� (x,y)
        parent_centre_pos   % ��Ҫ�洢umi����������
        sector_centre       % ��Ҫ�洢uma����������
        femto_room          % �洢femto�з����С
        O_InH;
        sectors             % eNodeB�ϵ�С��ʵ��
        servering_height    % ��վ����߶�
        in_interf_eNodeB_sites    % ��һ��ϵͳ���Ż�վ
        ad_interf_eNodeB_sites       % �ڶ���ϵͳ���Ż�վ
        name
        altitude = 0;
        site_name
        
        site_type           
        clock               % ����ʱ��
    end

    methods
        function print(obj)
            fprintf('eNodeB %d, position (x,y): (%3.2f,%3.2f), altitude: %3.0fm, %d attached UEs\n',obj.id,obj.pos(1),obj.pos(2),obj.altitude,obj.attached_UEs);
            fprintf('Sector: ');
            
            fprintf('\n');
            for s_=1:length(obj.sectors)
                obj.sectors(s_).print;
            end
        end
        
      
        function clear(obj)
            obj.sectors          = [];
            obj.in_interf_eNodeB_sites = [];
            obj.ad_interf_eNodeB_sites = [];
            obj.clock            = [];
        end
        
         %% �ж��Ƿ����û����ӵ���eNodeB
        function is_attached = userIsAttached(obj,user)
            for s_ = 1:length(obj.sectors)
                if obj.sectors(1).userIsAttached(user);
                    is_attached = true;
                    return
                end
            end
            is_attached = false;
        end
        %% �õ���eNodeB���ӵ��û��� 
        function number_of_atached_UEs = attached_UEs(obj)
            temp = zeros(1,length(obj.sectors));
            for s_ = 1:length(obj.sectors)
                temp(s_) = obj.sectors(s_).attached_UEs;
            end
            number_of_atached_UEs = sum(temp);
        end
    end
end