function NR_calculate_cell_coverage(SYS_config,networkPathlossMap,eNodeB_sites,eNodeB_sectors,networkShadowFadingMap)
%����ϵͳ���ó�������С������
%���������
%SYS_config��LTE��������ϵͳ
%networkPathlossMap��·��ͼ��
%eNodeB_sites����վ
%eNodeBs_sectors������
%networkShadowFadingMap����Ӱ˥��ӳ�����ͼ

if SYS_config.debug_level>=1
    fprintf('Calculating cell capacity:');
end

%% ��SNR
num_eNodeBs = length(eNodeB_sites);

shift_pathloss = networkPathlossMap.pathloss;
penetration_loss = networkPathlossMap.penetration_loss;
num_hetnet_sites = networkPathlossMap.num_hetnet_sites;
num_hetnet_sectors = networkPathlossMap.num_hetnet_sectors;
% UE�����������棬BS��û���������棬���ǽ���ʱ��Ȼ����BS��������
RX_powers_W = zeros(size(networkPathlossMap.pathloss));% ��ʼ�����չ��ʾ���
for b_ = 1:num_eNodeBs
    
    % ��Ӱ˥��
    shadow_fading_current_eNodeB = 10.^(networkShadowFadingMap.pathloss(:,:,b_)/10);
    
    for s_ = 1:length(eNodeB_sites(b_).sectors)
        s_idx = eNodeB_sites(b_).sectors(s_).eNodeB_id;
        %������չ���
        RX_powers_W(:,:,s_idx) = eNodeB_sites(b_).sectors(s_).max_power./10.^(shift_pathloss(:,:,s_idx)/10) ./ shadow_fading_current_eNodeB ./ 10.^(penetration_loss(:,:,s_idx)/10);
    end
end
%% ����SNR
% ���ﲻ��һ����վͬʱ������ٸ�UE������Ϊ��վ��һ��UEʹ��ȫ����Դ����Ϊ��������Ӱ�����վ�ķ�������
if SYS_config.debug_level>=1
    fprintf('calculating SNR-->');
end
SNR_linear_all = zeros(size(RX_powers_W));
thermal_noise_W = 10^(SYS_config.UE.thermal_noise_density/10) / 1000 * SYS_config.bandwidth * 10^(SYS_config.UE_receiver_noise_figure/10);

tot_eNodeBs = size(RX_powers_W,3);%����άΪС������
if SYS_config.isDouble
    num_first_sectors = networkPathlossMap.num_first_sectors;
    num_first_sites = networkPathlossMap.num_first_sites;
end
for s_=1:tot_eNodeBs
    SNR_linear_all(:,:,s_) = RX_powers_W(:,:,s_)/thermal_noise_W;
end
SNR_dB_all = 10*log10(SNR_linear_all);

% ��ÿ�����ص���յ���SNR��������ѡ����յ����SNRʱ��Ӧ�Ļ�վ��Ϊ�����վ���н���
if SYS_config.isDouble
    [~,SNR_IX] = sort(SNR_dB_all(:,:,1:num_first_sectors),3);
    [~,SNR_IX2] = sort(SNR_dB_all(:,:,num_first_sectors+1:tot_eNodeBs),3);% +num_first_sectors����Ϊ��������������Ǵ�1��ʼ��
else
    [~,SNR_IX] = sort(SNR_dB_all,3);
end

SNR_assignment  = SNR_IX(:,:,end);%SNR_assignment����ÿ�����ص�����SNR��Ӧ���������
if SYS_config.isDouble
    SNR_assignment2 = SNR_IX2(:,:,end)+num_first_sectors;% +num_first_sectors����Ϊ��������������Ǵ�1��ʼ��
else
    SNR_assignment2 = [];
end

% 3dB hand over
if SYS_config.hand_over
    if SYS_config.debug_level>=1
        fprintf('implement 3dB handover-->');
    end
    if SYS_config.isDouble
        % ��һ��ϵͳ
        SNR_dB_all_tmp = SNR_dB_all(:,:,1:num_first_sectors);
        for s_ = 1:num_first_sectors
            [m,n] = find(SNR_assignment == s_);% ��ǰ��վ��������ص�
            SNR_dB_all_tmp(m,n,s_) = SNR_dB_all_tmp(m,n,s_) - 3;% ��Щ���ص�SNR-3
            [~,SNR_IX_tmp] = sort(SNR_dB_all_tmp,3);% -3����������
            
            for i = 1:length(m)
                if SNR_IX_tmp(m(i),n(i),end) ~= s_ % ����Щ���ص���յ����SNR�Ļ�վ�����˸ı�
                    index = find(SNR_dB_all_tmp(m(i),n(i),:)>=SNR_dB_all_tmp(m(i),n(i),s_));% �ҳ���-3��SNRֵ���ȵĻ�վ
                    select = index(NR_randperm(length(index),1));% ����Щ��վ�����ѡ��һ����Ϊ�����վ
                    SNR_assignment(m(i),n(i)) = select;
                end
            end
        end
        % �ڶ���ϵͳ
        SNR_dB_all_tmp = SNR_dB_all(:,:,num_first_sectors+1:tot_eNodeBs);
        for s_ = 1:tot_eNodeBs-num_first_sectors
            [m,n] = find(SNR_assignment2 == (s_+num_first_sectors));
            SNR_dB_all_tmp(m,n,s_) = SNR_dB_all_tmp(m,n,s_) - 3;
            [~,SNR_IX2_tmp] = sort(SNR_dB_all_tmp,3);
            
            for i = 1:length(m)
                if SNR_IX2_tmp(m(i),n(i),end) ~= s_
                    index = find(SNR_dB_all_tmp(m(i),n(i),:)>=SNR_dB_all_tmp(m(i),n(i),s_));
                    select = index(NR_randperm(length(index),1));
                    SNR_assignment2(m(i),n(i)) = select+num_first_sectors;
                end
            end
        end
        
    else
        SNR_dB_all_tmp = SNR_dB_all;
        for s_ = 1:tot_eNodeBs
            [m,n] = find(SNR_assignment == s_);% ��ǰ��վ��������ص�
            SNR_dB_all_tmp(m,n,s_) = SNR_dB_all_tmp(m,n,s_) - 3;% ��Щ���ص�SNR-3
            [~,SNR_IX_tmp] = sort(SNR_dB_all_tmp,3);% -3����������
            
            for i = 1:length(m)
                if SNR_IX_tmp(m(i),n(i),end) ~= s_ % ����Щ���ص���յ����SNR�Ļ�վ�����˸ı�
                    index = find(SNR_dB_all_tmp(m(i),n(i),:)>=SNR_dB_all_tmp(m(i),n(i),s_));% �ҳ���-3��SNRֵ���ȵĻ�վ
                    select = index(NR_randperm(length(index),1));% ����Щ��վ�����ѡ��һ����Ϊ�����վ
                    SNR_assignment(m(i),n(i)) = select;
                end
            end
        end
        SNR_assignment2 = [];
    end
   
end

% ����õ���SNR_assignment��˼������SNR_assignment=(1,1,2,3,3,3),����1,1������ط��ɻ�վ1���񣬣�1,5������ط��ɻ�վ3����

% Ϊ�˽�UE�����᰸�ϵ����˷�Χ�ڣ��޸�SINR_assignment
% ����˼·�ǣ��Ƚ�ԭ����SINR_assignment����ú�-1��-1�������ڸô�����
% Ȼ���ҳ���������ĵط�����Ӧ����0
% �����SINR_assignment����0�ĵط���ԭ����ֵ�滻
if SYS_config.debug_level>=1
    fprintf('restrict cell coverage\n');
end
switch SYS_config.scene_type
    case {'UMA','RMa'}
        ISD = SYS_config.ISD;
        SNR_assignment_tmp = SNR_assignment;% �ȱ���ԭ����SNR_assignment
        SNR_assignment(:,:) = -1;% ��SNR_assignment����-1��-1�������ڴ˴��������˼
        distance_matrix = networkPathlossMap.distance_matrix;
        if SYS_config.isDouble
            n_enb = num_eNodeBs/2;
        else
            n_enb = num_eNodeBs;
        end
        for b_ = 1:n_enb
            [m,n] = find(distance_matrix(:,:,b_)<=2*ISD/3 & distance_matrix(:,:,b_)>=35);%UMA�����λ�ò����ڻ�վ��35m��Χ��
            for i = 1:length(m)
                for s_ = 1:length(eNodeB_sites(b_).sectors)
                    isinhex = isInHex([n(i),m(i)],eNodeB_sites(b_).sector_centre(s_).pos,ISD, networkPathlossMap);
                    if isinhex
                       SNR_assignment(m(i),n(i)) = 0; % �ܹ�����ĵط���0
                    end
                end
            end
        end
        [m,n] = find(SNR_assignment(:,:)==0);
        for i = 1:length(m)
            SNR_assignment(m(i),n(i)) = SNR_assignment_tmp(m(i),n(i));% ����������ĵط���ԭ����ֵ�滻
        end
        if SYS_config.isDouble
            SNR_assignment_tmp2 = SNR_assignment2;
            SNR_assignment2(:,:) = -1;
            for b_ = n_enb+1:num_eNodeBs
                [m,n] = find(distance_matrix(:,:,b_)<=2*ISD/3 & distance_matrix(:,:,b_)>=35);
                for i = 1:length(m)
                    for s_ = 1:length(eNodeB_sites(b_).sectors)
                        isinhex = isInHex([n(i),m(i)],eNodeB_sites(b_).sector_centre(s_).pos,ISD, networkPathlossMap);
                        if isinhex
                            SNR_assignment2(m(i),n(i)) = 0;
                        end
                    end
                end
            end
            [m,n] = find(SNR_assignment2(:,:)==0);
            for i = 1:length(m)
                SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
            end
        end
    case 'UMI'
        if ~SYS_config.isManhattan
            UMi_r = SYS_config.UMi_r;
            SNR_assignment_tmp = SNR_assignment;
            SNR_assignment(:,:) = -1;
            distance_matrix = networkPathlossMap.distance_matrix;
            micro_centre_distance_matrix = networkPathlossMap.micro_centre_distance_matrix;
            if SYS_config.isDouble
                n_enb = num_eNodeBs/2;
            else
                n_enb = num_eNodeBs;
            end
            for b_ = 1:n_enb
                [m,n] = find(micro_centre_distance_matrix(:,:,b_)<=UMi_r & distance_matrix(:,:,b_)>=3);%UMI���������28.9��ΧȦ�ڣ��Ҳ�����΢վ��3m��Χ�� 
                for i = 1:length(m)
                        SNR_assignment(m(i),n(i)) = 0;                    
                end
            end
            [m,n] = find(SNR_assignment(:,:)==0);
            for i = 1:length(m)
                SNR_assignment(m(i),n(i)) = SNR_assignment_tmp(m(i),n(i));
            end
            if SYS_config.isDouble
                SNR_assignment_tmp2 = SNR_assignment2;
                SNR_assignment2(:,:) = -1;
                for b_ = n_enb+1:num_eNodeBs
                    [m,n] = find(micro_centre_distance_matrix(:,:,b_)<=UMi_r & distance_matrix(:,:,b_)>=3);
                    for i = 1:length(m)
                        SNR_assignment2(m(i),n(i)) = 0;%������ŵ��ڻ�վ���
                    end
                end
                [m,n] = find(SNR_assignment2(:,:)==0);
                for i = 1:length(m)
                    SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
                end
            end
        else
            % do nothing
        end
    case {'InH','InH2'}
        % do nothing
    case 'InH3'
        ISD = SYS_config.ISD;
        SNR_assignment_tmp = SNR_assignment;% �ȱ���ԭ����SNR_assignment
        SNR_assignment(:,:) = -1;% ��SNR_assignment����-1��-1�������ڴ˴��������˼
        distance_matrix = networkPathlossMap.distance_matrix;
        if SYS_config.isDouble
            n_enb = num_eNodeBs/2;
        else
            n_enb = num_eNodeBs;
        end
        
        [m,n] = size(SNR_assignment);%UMA�����λ�ò����ڻ�վ��35m��Χ��
        for i = 1:m
            for j = 1:n
                pos=NR_common_pixel_to_pos([i j],[networkPathlossMap.roi_y(1) networkPathlossMap.roi_x(1)],networkPathlossMap.data_res);
                if pos(2)>=0 && pos(2)<=120 && pos(1)>=0 && pos(1)<=50
                    inarea = 1;
                else
                    inarea = 0;
                end
                if inarea
                    SNR_assignment(i,j) = 0; % �ܹ�����ĵط���0
                end
            end
        end
        
        [m,n] = find(SNR_assignment(:,:)==0);
        for i = 1:length(m)
            SNR_assignment(m(i),n(i)) = SNR_assignment_tmp(m(i),n(i));% ����������ĵط���ԭ����ֵ�滻
        end
        if SYS_config.isDouble
            SNR_assignment_tmp2 = SNR_assignment2;
            SNR_assignment2(:,:) = -1;
            [m,n] = size(SNR_assignment2);
            for i = 1:m
                for j = 1:n
                    pos=NR_common_pixel_to_pos([i j],[networkPathlossMap.roi_y(1) networkPathlossMap.roi_x(1)],networkPathlossMap.data_res);
                    if pos(2)>=0 && pos(2)<=120 && pos(1)>=50+SYS_config.InH3_d && pos(1)<=100+SYS_config.InH3_d
                        inarea = 1;
                    else
                        inarea = 0;
                    end
                    if inarea
                        SNR_assignment2(i,j) = 0; % �ܹ�����ĵط���0
                    end
                end
            end
            [m,n] = find(SNR_assignment2(:,:)==0);
            for i = 1:length(m)
                SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
            end
        end
    case 'UMa_to_UMi'
        if ~SYS_config.isManhattan
            % ΢С����Ϊ������
            UMi_r = SYS_config.UMi_r;
            SNR_assignment_tmp = SNR_assignment;
            SNR_assignment(:,:) = -1;
            distance_matrix = networkPathlossMap.distance_matrix;
            micro_centre_distance_matrix = networkPathlossMap.micro_centre_distance_matrix;
            for b_ = 1:num_first_sites
                [m,n] = find(micro_centre_distance_matrix(:,:,b_)<=UMi_r & distance_matrix(:,:,b_)>=3);
                for i = 1:length(m)
                    SNR_assignment(m(i),n(i)) = 0;
                end
            end
            [m,n] = find(SNR_assignment(:,:)==0);
            for i = 1:length(m)
                SNR_assignment(m(i),n(i)) = SNR_assignment_tmp(m(i),n(i));
            end
        else
            % ΢С��Ϊ������
            distance_matrix = networkPathlossMap.distance_matrix;
            O_M = [-540,-540];
            K_M = [540,540];
            O_M_pix = NR_common_pos_to_pixel(O_M, networkPathlossMap.coordinate_origin, networkPathlossMap.data_res);
            K_M_pix = NR_common_pos_to_pixel(K_M, networkPathlossMap.coordinate_origin, networkPathlossMap.data_res);
            M_row = O_M_pix(2):K_M_pix(2);
            M_col = O_M_pix(1):K_M_pix(1);
            tmp = zeros(size(SNR_assignment));
            tmp(M_row,M_col) = 1;
            SNR_assignment = SNR_assignment .* tmp;
        end
        
        ISD = SYS_config.ISD;
        SNR_assignment_tmp2 = SNR_assignment2;
        SNR_assignment2(:,:) = -1;
        for b_ = num_first_sites+1:num_eNodeBs
            [m,n] = find(distance_matrix(:,:,b_)<=2*ISD/3 & distance_matrix(:,:,b_)>=35);
            for i = 1:length(m)
                for s_ = 1:length(eNodeB_sites(b_).sectors)
                    isinhex = isInHex([n(i),m(i)],eNodeB_sites(b_).sector_centre(s_).pos,ISD, networkPathlossMap);
                    if isinhex
                        SNR_assignment2(m(i),n(i)) = 0;
                    end
                end
            end
        end
        [m,n] = find(SNR_assignment2(:,:)==0);
        for i = 1:length(m)
            SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
        end
    case {'UMa_to_InH','UMi_to_InH'}
        % ��O_InH
        O_InH = eNodeB_sites(1).O_InH;
        if ~SYS_config.isFemto
            K_InH = [O_InH(1)+120,O_InH(2)+50];
        else
            K_InH = [O_InH(1)+120,O_InH(2)+40];
        end
        O_InH_pix = NR_common_pos_to_pixel(O_InH, networkPathlossMap.coordinate_origin, networkPathlossMap.data_res);
        K_InH_pix = NR_common_pos_to_pixel(K_InH, networkPathlossMap.coordinate_origin, networkPathlossMap.data_res);
        InH_row = O_InH_pix(2):K_InH_pix(2);
        InH_col = O_InH_pix(1):K_InH_pix(1);
        % ��inh SINR_assignment
        tmp = zeros(size(SNR_assignment));
        tmp(InH_row,InH_col) = 1;
        SNR_assignment = SNR_assignment .* tmp;
        
        switch SYS_config.scene_type
            case 'UMa_to_InH'
                ISD = SYS_config.ISD;
                SNR_assignment_tmp2 = SNR_assignment2;
                SNR_assignment2(:,:) = -1;
                distance_matrix = networkPathlossMap.distance_matrix;
                for b_ = num_first_sites+1:num_eNodeBs
                    [m,n] = find(distance_matrix(:,:,b_)<=2*ISD/3 & distance_matrix(:,:,b_)>=35);
                    for i = 1:length(m)
                        for s_ = 1:length(eNodeB_sites(b_).sectors)
                            isinhex = isInHex([n(i),m(i)],eNodeB_sites(b_).sector_centre(s_).pos,ISD, networkPathlossMap);
                            if isinhex
                                SNR_assignment2(m(i),n(i)) = 0;
                            end
                        end
                    end
                end
                [m,n] = find(SNR_assignment2(:,:)==0);
                for i = 1:length(m)
                    SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
                end
            case 'UMi_to_InH'
                UMi_r = SYS_config.UMi_r;
                SNR_assignment_tmp2 = SNR_assignment2;
                SNR_assignment2(:,:) = -1;
                distance_matrix = networkPathlossMap.distance_matrix;
                micro_centre_distance_matrix = networkPathlossMap.micro_centre_distance_matrix;
                for b_ = num_first_sites+1:num_eNodeBs
                    [m,n] = find(micro_centre_distance_matrix(:,:,b_)<=UMi_r & distance_matrix(:,:,b_)>=3);
                    for i = 1:length(m)
                        SNR_assignment2(m(i),n(i)) = 0;
                    end
                end
                [m,n] = find(SNR_assignment2(:,:)==0);
                for i = 1:length(m)
                    SNR_assignment2(m(i),n(i)) = SNR_assignment_tmp2(m(i),n(i));
                end
        end
end

% ����������С
if SYS_config.isDouble
    n_sector1 = num_first_sectors;
    n_sector2 = num_hetnet_sectors;
else
    n_sector1 = length(eNodeB_sectors);
end
cell_sizes = zeros(1,n_sector1);
cell_centers_pixel = zeros(n_sector1,2);

for s_idx = 1:n_sector1
    cell_sizes(s_idx) = sum(SNR_assignment(:)==s_idx);% ����SNR_assignment=(1,1,2,3,3,3),sum(SNR_assignment(:)==1)�õ�2����˼�ǻ�վ1�������������ص�
    [row,col] = find(SNR_assignment==s_idx);
    cell_centers_pixel(s_idx,:) = [mean(col) mean(row)];
end
if SYS_config.isDouble% �ڶ���ϵͳ
    cell_sizes2 = zeros(1,n_sector2);
    cell_centers_pixel2 = zeros(n_sector2,2);
    for s_idx = n_sector1+1:length(eNodeB_sectors)
        cell_sizes2(s_idx-n_sector1) = sum(SNR_assignment2(:)==s_idx);
        [row,col] = find(SNR_assignment2==s_idx);
        cell_centers_pixel2(s_idx-n_sector1,:) = [mean(col) mean(row)];
    end
    cell_centers2 = NR_common_pixel_to_pos(cell_centers_pixel2,networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);% �ڶ���ϵͳ����������
else
    cell_sizes2 = [];
    cell_centers2 = [];
end
% ��һ��ϵͳ����������
cell_centers = NR_common_pixel_to_pos(cell_centers_pixel,networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);% ��¼�������ģ��ڻ�ͼʱ��ǳ��������

%% ��������
networkPathlossMap.sector_assignment = SNR_assignment;
networkPathlossMap.sector_assignment_double = SNR_assignment2;
networkPathlossMap.sector_sizes = cell_sizes;% ��������������ص�����һ��1*n_sector1ά����
networkPathlossMap.sector_sizes_double = cell_sizes2;
networkPathlossMap.sector_centers = cell_centers;
networkPathlossMap.sector_centers_double = cell_centers2;

end

% �ж�ĳ���Ƿ�����������
function isinhex = isInHex( dot_pix,sector_centre,ISD, networkPathlossMap)

dot_pos = NR_common_pixel_to_pos(dot_pix,networkPathlossMap.coordinate_origin,networkPathlossMap.data_res);
rot = @(w_) [cos(w_),-sin(w_);sin(w_),cos(w_)];
for i_ = 1:7
    tmp_hex(i_,:) = sector_centre + (ISD/3*(rot(pi/3)^(i_-1)*[1;0]).'); % ��[1;0]��ʾֻȡ��һ��
end
isinhex = inpolygon(dot_pos(1),dot_pos(2),tmp_hex(:,1),tmp_hex(:,2));

end
