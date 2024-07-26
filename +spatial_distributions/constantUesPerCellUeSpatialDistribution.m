classdef constantUesPerCellUeSpatialDistribution < spatial_distributions.generalUeSpatialDistribution
    % ���������һ����

    properties
        % ÿ��С��UE��
        UE_per_eNodeB
        % �칹�µڶ���ϵͳÿ��С��UE��
        UE_per_eNodeB2
        % ˫ϵͳʱ��һ��ϵͳ������
        num_first_sectors
    end
    
    methods
        % Class constructor.
        function obj = constantUesPerCellUeSpatialDistribution(networkPathlossMap,eNodeB_sites,eNodeB_sectors,SYS_config)
            % Fill in basic parameters (handled by the superclass constructor)
            obj = obj@spatial_distributions.generalUeSpatialDistribution(networkPathlossMap,eNodeB_sites,eNodeB_sectors);%@����������Ĺ��캯����@���ڵ��û����Ա����
            switch SYS_config.scene_type
                case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
                    obj.UE_per_eNodeB = SYS_config.UE_per_eNodeB;
                    obj.UE_per_eNodeB2 = SYS_config.UE_per_eNodeB2;
                    obj.num_first_sectors = length(eNodeB_sectors)-networkPathlossMap.num_hetnet_sectors;% ֻ����˫ϵͳ��
                otherwise
                    obj.UE_per_eNodeB = SYS_config.UE_per_eNodeB;
                    obj.UE_per_eNodeB2 = 0;
                    obj.num_first_sectors = length(eNodeB_sectors)/2;% ֻ����˫ϵͳ��
            end
        end
        
        function user_positions_pixels = generate_UE_positions(obj,SYS_config)
            networkPathlossMap = obj.networkPathlossMap;
            sector_surfaces    = networkPathlossMap.sector_sizes;
            sector_surfaces2 = networkPathlossMap.sector_sizes_double;
            eNodeB_sectors            = obj.eNodeB_sectors;%sector
            
            % ���ֻ�վ�ķ���ΧΪ0������������������һ�ֿ����Ƿֱ���̫�֣����п�����ǰ�����﷢���˴���
            % �������������ʱ�ᵼ����SINR���ﱨ��
            if sum(sector_surfaces(:)==0)>0
                warning('Some sector_surfaces are zero. No UEs will be generated there. Maybe a too big map resolution value?');
            end
            
            if SYS_config.isDouble
                if sum(sector_surfaces2(:)==0)>0
                    warning('Some sector_surfaces2 are zero. No UEs will be generated there. Maybe a too big map resolution value?');
                end
            end
            
            % users_sector��ĳ������Ҫ����UE��
            norm_sector_surface = ones(size(sector_surfaces));
            norm_sector_surface2 = ones(size(sector_surfaces2));
            switch SYS_config.scene_type
                case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}%ûŪ��
                    users_sector = round(norm_sector_surface*obj.UE_per_eNodeB).*(sector_surfaces>0);
                    users_sector2 = round(norm_sector_surface2*obj.UE_per_eNodeB2).*(sector_surfaces2>0);
                otherwise
                    users_sector = round(norm_sector_surface*obj.UE_per_eNodeB).*(sector_surfaces>0);
                    users_sector2 = round(norm_sector_surface2*obj.UE_per_eNodeB).*(sector_surfaces2>0);
            end
            
            % sector_positions������������������ص��λ�ã�user_positions_pixels_per_cell��������������ص��λ��
            sector_positions = cell(size(networkPathlossMap.sector_sizes));%cell��һ�ֿ��Դ洢��ͬ�������ݵ�����
            user_positions_pixels_per_cell = cell(size(networkPathlossMap.sector_sizes));
            sector_positions2 = cell(size(networkPathlossMap.sector_sizes_double));
            user_positions_pixels_per_cell2 = cell(size(networkPathlossMap.sector_sizes_double));
            
            % Assign random positions to each UE
            if SYS_config.isDouble
                n_sector = obj.num_first_sectors;
            else
                n_sector = length(eNodeB_sectors);
            end
            for s_idx = 1:n_sector
                if users_sector(s_idx)~=0                   
                    [row,col] = find(networkPathlossMap.sector_assignment(:,:)==s_idx);
                    sector_positions{s_idx} = [col,row];%cell�Ĵ����ű�ʾ����
                    % �ڻ�վ�ķ���Χ�����ѡ��UE_per_eNodeB������ΪUE��λ��
                    user_positions_pixels_per_cell{s_idx} = sector_positions{s_idx}(randi(size( sector_positions{s_idx},1)  ,  [1 users_sector(s_idx)]  ),:);
                else
                    sector_positions{s_idx}      = [];
                    user_positions_pixels_per_cell{s_idx} = [];
                end
            end
            if SYS_config.isDouble % �ڶ���ϵͳ���ò�ͬ��sector_assignment
                for s_idx = n_sector+1:length(eNodeB_sectors)
                    if users_sector2(s_idx-n_sector)~=0
                        [row,col] = find(networkPathlossMap.sector_assignment_double(:,:)==s_idx);
                        sector_positions2{s_idx-n_sector} = [col,row];%cell�Ĵ����ű�ʾ����
                        % �ڻ�վ�ķ���Χ�����ѡ��UE_per_eNodeB������ΪUE��λ��
                        user_positions_pixels_per_cell2{s_idx-n_sector} = sector_positions2{s_idx-n_sector}(randi(size( sector_positions2{s_idx-n_sector},1)  ,  [1 users_sector2(s_idx-n_sector)]  ),:);
                    else
                        sector_positions2{s_idx-n_sector}      = [];
                        user_positions_pixels_per_cell2{s_idx-n_sector} = [];
                    end
                end
                user_positions_pixels_per_cell = [user_positions_pixels_per_cell user_positions_pixels_per_cell2];
            end
            
            % �����ⲽ�ǽ�����õ���UE��λ�����ϳ�n_UE�У�2�е���ʽ����һ����UE��x�ᣬ�ڶ�����UE��y��
            user_positions_pixels_total_UEs = 0;
            for i_=1:length(user_positions_pixels_per_cell)
                if ~isempty(user_positions_pixels_per_cell{i_})
                    user_positions_pixels_total_UEs = user_positions_pixels_total_UEs + size(user_positions_pixels_per_cell{i_},1);
                end
            end
            user_positions_pixels = zeros(user_positions_pixels_total_UEs,2);
            l_=1;
            for i_=1:length(user_positions_pixels_per_cell)
                for j_=1:size(user_positions_pixels_per_cell{i_},1)
                    user_positions_pixels(l_,:) = user_positions_pixels_per_cell{i_}(j_,:);
                    l_=l_+1;
                end
            end
        end
    end
    
end

