classdef generalUeSpatialDistribution < handle
    % �����ǿ���UE����ĸ���
    % ���캯�������ǰ�networkPathlossMap��eNodeB_sites��eNodeB_sectors����ȥ
    % �������ǰ��ƽ̨�̳���������ǰ��ƽ̨���б�����㷽ʽ�������ǵ�ƽ̨��ֻ����ÿ��С������ͬ��UE������
    
    properties
        networkPathlossMap % ��Ҫ�ǵ�����sector_assignment����
        eNodeB_sites              % ��վ
        eNodeB_sectors            % ����
    end
    
    methods
        function obj = generalUeSpatialDistribution(networkPathlossMap,eNodeB_sites,eNodeB_sectors)
            obj.networkPathlossMap = networkPathlossMap;
            obj.eNodeB_sites              = eNodeB_sites;
            obj.eNodeB_sectors            = eNodeB_sectors;
        end
    end

    methods(Abstract)
        user_positions_pixels = generate_UE_positions(obj,varargin)
    end
    
end

