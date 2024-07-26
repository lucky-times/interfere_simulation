function SYS_config = NR_load_params_again( SYS_config )
% ��ֹĳЩ�������ó�����Щ����Ӱ�쵽���˵Ľ����Ͳ��������ļ�������
% ���������
% SYS_config����������
%
% ���������
% SYS_config����������

% check
try
    if SYS_config.n_snapshot~=SYS_config.n_snapshot2 % �ݲ�֧���칹ʱ����С����ͬUE��Ŀ�����
        SYS_config.n_snapshot2 = SYS_config.n_snapshot;
    end
catch
    % ��SYS_config.n_snapshot2û�ж��壬˵�������칹����
    % do nothing
end

% �Բ�ͬ�����Ĳ�����������
switch SYS_config.scene_type
    case 'UMA'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;% ���������վͬʱ�����UE����˵õ�һ��С����UE��
        SYS_config.isManhattan = false;
        SYS_config.isFemto = false;
        SYS_config.isS2F = true;
    case 'UMI'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.shift_mode = 1;
        SYS_config.isNewUMA = false;
        SYS_config.isWraparound = false;
        SYS_config.isFemto = false;
        SYS_config.isS2F = true;
    case {'InH','InH2','InH3'}
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.ISD = [];
        SYS_config.shift_mode = 0;
        SYS_config.isNewUMA = false;
        SYS_config.isManhattan = false;
        SYS_config.isWraparound = false;
        if SYS_config.isFemto
            SYS_config.shift_mode = 1;
        end
        SYS_config.isS2F = true;
    case 'RMa'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.isManhattan = false;
        SYS_config.isFemto = false;
        SYS_config.isS2F = true;
    case 'UMa_to_UMi'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.UE_per_eNodeB2 = SYS_config.n_snapshot2*SYS_config.n_UE_served_per_BS;
        SYS_config.isDouble = true;
        SYS_config.isNewUMA = true;
        SYS_config.isWraparound = false;
        SYS_config.isFemto = false;
        if SYS_config.isManhattan
            SYS_config.ISD = 500;
            SYS_config.map_resolution = 15;
            SYS_config.shadow_fading_map_resolution=15;
        end
    case 'UMa_to_InH'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.UE_per_eNodeB2 = SYS_config.n_snapshot2*SYS_config.n_UE_served_per_BS;
        SYS_config.isDouble = true;
        SYS_config.isNewUMA = false;
        SYS_config.isManhattan = false;
        SYS_config.isWraparound = false;
    case 'UMi_to_InH'
        SYS_config.UE_per_eNodeB = SYS_config.n_snapshot*SYS_config.n_UE_served_per_BS;
        SYS_config.UE_per_eNodeB2 = SYS_config.n_snapshot2*SYS_config.n_UE_served_per_BS;
        SYS_config.isDouble = true;
        SYS_config.isNewUMA = false;
        SYS_config.isManhattan = false;
        SYS_config.isWraparound = false;
end

% ������ļ�����
name_generator = utils.naming;
SYS_config.results_file = name_generator.results_file(SYS_config);
SYS_config.beam_cache_file = name_generator.beam_cache_file(SYS_config);

end

