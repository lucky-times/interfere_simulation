classdef generalPathlossModel < handle
% ����·��ģ�������ߵ�Ԫ����
    properties
        name
    end
    methods (Abstract)
        pathloss_in_db = pathloss(distance)
    end
    methods (Static)
        function calculate_pathloss_maps(SYS_config,eNodeB_sites,networkMacroscopicPathlossMap)
            % ��һЩ������ʼ��
            total_sectors         = length([eNodeB_sites.sectors]);
            data_res              = networkMacroscopicPathlossMap.data_res;
            roi_x                 = networkMacroscopicPathlossMap.roi_x;
            roi_y                 = networkMacroscopicPathlossMap.roi_y;
            roi_maximum_pixels    = NR_common_pos_to_pixel( [roi_x(2) roi_y(2)], [roi_x(1) roi_y(1)], data_res);
            roi_height_pixels     = roi_maximum_pixels(2);
            roi_width_pixels      = roi_maximum_pixels(1);
            distance_matrix       = zeros(roi_height_pixels,roi_width_pixels,length(eNodeB_sites));
            cell_pathloss_data    = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            sector_antenna_gain   = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            % ����ά��ͼ�йر���
            UE_antenna_gain_1     = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            UE_antenna_gain_2     = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            UE_antenna_gain     = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            penetration_loss = zeros(roi_height_pixels,roi_width_pixels,total_sectors);
            d_2d_in = zeros(roi_height_pixels,roi_width_pixels);
            
            site_positions     = reshape([eNodeB_sites.pos],2,[])';
            site_positions_pix = NR_common_pos_to_pixel( site_positions, [roi_x(1) roi_y(1)], data_res);
            %% ����3ά��ͼ
            switch SYS_config.scene_type
                case {'UMA' , 'UMI','UMa_to_UMi'}
                    [row,col] = size(d_2d_in);
                    if ~SYS_config.isManhattan %UMA��UMi��Rand dropģ��
                        % 'UMI','UMa_to_UMi'��΢С������������ִ����һ��
%                         if ~SYS_config.isNewUMA %��ǰ��UMA������500ISD�������������Ϊ20%��20%
%                             index = NR_randperm(row*col,round(row*col*0.2));% 20%������
%                         else
%                             index = NR_randperm(row*col,round(row*col*0));% 0%������
%                         end
                        switch SYS_config.scene_type
                            case 'UMA'
                                index = NR_randperm(row*col,round(row*col*0));
                            case 'UMI'
                                index = NR_randperm(row*col,round(row*col*0.2));
                            case 'UMa_to_UMi'
                                index = NR_randperm(row*col,round(row*col*0.2));
                        end
                        for i = 1:length(index)
                            if index == 0
                                d_2d_in = zeros(size(d_2d_in));
                            else
                                d_2d_in(index(i)) = min(25*rand(1,1),25*rand(1,1)); %����d_2d_in����������ľ��룩�����������û���������·�����
                            end
                        end
                    else % ������
                        % 'UMI','UMa_to_UMi'��΢С����������ִ����һ��
                        switch SYS_config.scene_type
                            case 'UMI'
                                % 'UMI'��΢С����������
                                % ����ͼ�ֱ���Ϊ5����һ�����ص�ռ�ľ���Ϊ5m
                                % ��Ϊ¥��ͽֵ���������ԣ�SYS_config.map_resolutionӦ����һ��������75��15����
                                % ����������������������ϵͳ�Լ��Ĺ���û��Ӱ�죬���ǶԺ�����������ٵ��칹�ᱨ����Ϊ�����ģ��һ��
                                n_ = fix(75/SYS_config.map_resolution);% ¥��
                                m_ = fix(15/SYS_config.map_resolution);% �ֵ���
                                % �����½ǿ�ʼ������ӽ�����
                                for row_=m_+1:m_+n_ % ��һ�н�����
                                    d_2d_in(row_,1:n_) = min(25*rand(1,n_),25*rand(1,n_));
                                    for i=1:11
                                        d_2d_in(row_,i*n_+i*m_+1:(i+1)*n_+i*m_) = min(25*rand(1,n_),25*rand(1,n_));%�����������ֵ��ѡ���С��
                                    end
                                end
                                for j=1:11% ����11�н�����
                                    for row_=m_+j*n_+j*m_+1:m_+(j+1)*n_+j*m_
                                        d_2d_in(row_,1:n_) = min(25*rand(1,n_),25*rand(1,n_));
                                        for i=1:11
                                            d_2d_in(row_,i*n_+i*m_+1:(i+1)*n_+i*m_) = min(25*rand(1,n_),25*rand(1,n_));
                                        end
                                    end
                                end
                            case 'UMa_to_UMi'
                                % 'UMa_to_UMi'��΢С����������
                                % UMa�����µ�d_2d_in
                                if ~SYS_config.isNewUMA %��ǰ��UMA������500ISD�������������Ϊ80%��20%
                                    index = NR_randperm(row*col,round(row*col*0.8));% 80%������
                                else
                                    index = NR_randperm(row*col,round(row*col*0.2));% 20%������
                                end
                                for i = 1:length(index)
                                    d_2d_in(index(i)) = min(25*rand(1,1),25*rand(1,1)); %����d_2d_in����������ľ��룩�����������û���������·�����
                                end
                                % �����ٳ����µ�d_2d_in
                                O_M = [-540,-540];
                                K_M = [540,540];
                                O_M_pix = NR_common_pos_to_pixel(O_M, [roi_x(1) roi_y(1)], data_res);
                                K_M_pix = NR_common_pos_to_pixel(K_M, [roi_x(1) roi_y(1)], data_res);
                                M_row = O_M_pix(2):K_M_pix(2);
                                M_col = O_M_pix(1):K_M_pix(1);
                                d_2d_in_M = zeros(size(d_2d_in(M_row,M_col)));
                                
                                n_ = round(75/SYS_config.map_resolution);
                                m_ = round(15/SYS_config.map_resolution);
                                % �������������ˣ�����8��*�����һ��������
                                % ��ӽ�����ķ�������һ��һ�����ص�����
                                % һ�н��������˼������ͼ���������������ص㣨*������һ�н�����
                                % ****  ****
                                % ****  ****
                                % 
                                % ****  ****
                                % ****  ****
                                % �����½ǿ�ʼ������ӽ�����
                                for row_=m_+1:m_+n_ % ��һ�н�����
                                    d_2d_in_M(row_,1:n_) = min(25*rand(1,n_),25*rand(1,n_));
                                    for i=1:11
                                        d_2d_in_M(row_,i*n_+i*m_+1:(i+1)*n_+i*m_) = min(25*rand(1,n_),25*rand(1,n_));%�����������ֵ��ѡ���С��
                                    end
                                end
                                for j=1:11% ����11�н�����
                                    for row_=m_+j*n_+j*m_+1:m_+(j+1)*n_+j*m_
                                        d_2d_in_M(row_,1:n_) = min(25*rand(1,n_),25*rand(1,n_));
                                        for i=1:11
                                            d_2d_in_M(row_,i*n_+i*m_+1:(i+1)*n_+i*m_) = min(25*rand(1,n_),25*rand(1,n_));
                                        end
                                    end
                                end
                                % �����UMa�����µ�d_2d_in��Ȼ�����������µ�d_2d_in_M�滻d_2d_in��Ӧ�Ĳ���
                                d_2d_in(M_row,M_col) = d_2d_in_M;
                        end
                    end
                    % ue�߶�
                    h_ut = ones(size(d_2d_in)) * SYS_config.UE_height;
                    
                    % ���𡪡�һ��ʹ���һ��ʹ���

                    % �ʹ��� ������38.900��Ĵ���ʽ���㣩
                    L_glass = 2+0.2*SYS_config.frequency/1e9;
                    L_concrete = 5+4*SYS_config.frequency/1e9;
                    N_low = 0+4.4*randn(row,col);
                    low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*d_2d_in+N_low;
                    index = NR_randperm(row*col,round(row*col*0.5));%��·�𡢸�·���ȡһ��
                    low(index) = 0;
                    % �ߴ���
                    L_IRRglass = 23+0.3*SYS_config.frequency/1e9;
                    N_high = 0+6.5*randn(row,col);
                    high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*d_2d_in+N_high;
                    [m,n]=find(low(:,:)~=0);
                    for i = 1:length(m)
                        high(m(i),n(i)) = 0;
                    end
                    penetration = low+high;%�ܵĴ���
                    [m,n]=find(d_2d_in(:,:)==0);
                    for i = 1:length(m)
                        penetration(m(i),n(i)) = 0; %�����û�����Ϊ0
                    end
                    if SYS_config.isDouble% Ƶ�ʲ�һ��������һ��
                        % �ʹ���
                        L_glass = 2+0.2*SYS_config.frequency2/1e9;
                        L_concrete = 5+4*SYS_config.frequency2/1e9;
                        low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*d_2d_in+N_low;
                        low(index) = 0;
                        %�ߴ���
                        L_IRRglass = 23+0.3*SYS_config.frequency2/1e9;
                        high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*d_2d_in+N_high;
                        [m,n]=find(low(:,:)~=0);
                        for i = 1:length(m)
                            high(m(i),n(i)) = 0;
                        end
                        penetration2 = low+high;
                        [m,n]=find(d_2d_in(:,:)==0);
                        for i = 1:length(m)
                            penetration2(m(i),n(i)) = 0;
                        end
                    end
                case 'RMa'
                    % 50%���ڣ�50%���ϣ����ڵʹ��𣬳���û�д���ֻ��һ��N(��, ��P��
                    [row,col] = size(d_2d_in);
                    index = NR_randperm(row*col,round(row*col*0.5));% 50%������
                    for i = 1:length(index)
                        d_2d_in(index(i)) = min(25*rand(1,1),25*rand(1,1));
                    end
                    % ue�߶�
                    h_ut = ones(size(d_2d_in)) * 1.5;
                    
                    % �ʹ���                    
                    L_glass = 2+0.2*SYS_config.frequency/1e9;
                    L_concrete = 5+4*SYS_config.frequency/1e9;
                    N_low = 0+4.4*randn(row,col);
                    low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*d_2d_in+N_low;
                    penetration = low;
                    [m,n]=find(d_2d_in(:,:)==0);
                    for i = 1:length(m)
                        penetration(m(i),n(i)) = 9+5*randn;% ����û�д���ֻ��һ��N(��, ��P��
                    end
                    if SYS_config.isDouble
                        L_glass = 2+0.2*SYS_config.frequency2/1e9;
                        L_concrete = 5+4*SYS_config.frequency2/1e9;
                        low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*d_2d_in+N_low;
                        penetration2 = low;
                        [m,n]=find(d_2d_in(:,:)==0);
                        for i = 1:length(m)
                            penetration2(m(i),n(i)) = 9+5*randn;
                        end
                    end
            end
            
            % ��ʼ������ͽǶȾ���
            position_grid_pixels      = zeros(roi_height_pixels*roi_width_pixels,2);
            position_grid_pixels(:,1) = reshape(repmat(1:roi_width_pixels,roi_height_pixels,1),1,roi_width_pixels*roi_height_pixels);
            position_grid_pixels(:,2) = repmat(1:roi_height_pixels,1,roi_width_pixels)';
            position_grid_meters      = NR_common_pixel_to_pos(position_grid_pixels,networkMacroscopicPathlossMap.coordinate_origin,networkMacroscopicPathlossMap.data_res);
            
            % ��femto��qcol��qrow��qcol�Ǿ�����ֱ��ǽ�ĸ�����qrow�Ǿ���ˮƽ��ǽ�ĸ���
            if SYS_config.isFemto
                if ~isFemto_het
                    q_f_matrix = zeros(roi_height_pixels,roi_width_pixels,length(eNodeB_sites));
                    % q_col
                    % ÿ������Ĵ�С��20*20�����ǿ���ͨ�������㷨�жϸõ����ĸ�����
                    q_col_mod = mod(position_grid_meters(:,1),20);
                    index = find(q_col_mod == 0);
                    q_col = fix(position_grid_meters(:,1)/20)+1;
                    q_col(index) = q_col(index) - 1;
                    % q_row
                    q_row_mod = mod(position_grid_meters(:,2),20);
                    index = find(q_row_mod == 0);
                    q_row = fix(position_grid_meters(:,2)/20)+1;
                    q_row(index) = q_row(index) - 1;
                else
                    q_f_matrix = zeros(roi_height_pixels,roi_width_pixels,length(eNodeB_sites));
                    % q_col
                    position_col_matrix = reshape(position_grid_meters(:,1),roi_height_pixels,roi_width_pixels);
                    insite_col = position_col_matrix(InH_row,InH_col) + 60;% ��60��Ϊ�˵������꣬ʹ��ͷ��칹������һ�£������������
                    q_col_mod = mod(insite_col,20);
                    index = find(q_col_mod == 0);
                    q_col = fix(insite_col/20)+1;
                    q_col(index) = q_col(index) - 1;
                    % q_row
                    position_row_matrix = reshape(position_grid_meters(:,2),roi_height_pixels,roi_width_pixels);
                    insite_row = position_row_matrix(InH_row,InH_col) + 20;
                    q_row_mod = mod(insite_row,20);
                    index = find(q_row_mod == 0);
                    q_row = fix(insite_row/20)+1;
                    q_row(index) = q_row(index) - 1;
                end
            end
            
            %% ����С����·��
            s_idx = 1;
            all_sectors   = [eNodeB_sites.sectors];
            eNodeB_id_set = [all_sectors.eNodeB_id];
            % 2panel,+-180,UE��������
            rnd_phi_panel1=180-360*rand(size(h_ut)); %����-180~180�������ת��λ
            rnd_phi_panel2=rnd_phi_panel1-180;%����panel��180��
            for b_ = 1:length(eNodeB_sites)
                distances               = sqrt((position_grid_meters(:,1)-eNodeB_sites(b_).pos(1)).^2 + (position_grid_meters(:,2)-eNodeB_sites(b_).pos(2)).^2);%ÿ�����ص㵽enb�ľ���
                % ����ͼ��size�����������
                distance_matrix(:,:,b_) = reshape(distances,roi_height_pixels,roi_width_pixels);
                
                % �����ص�����վ���˶�����ǽ
                if SYS_config.isFemto
                    % �����վ���ĸ�����
                    if ~isFemto_het
                        site_pos = eNodeB_sites(b_).pos;
                    else
                        site_pos = [eNodeB_sites(b_).pos(1)+60,eNodeB_sites(b_).pos(2)+20];% �칹��������Ҫ��������ʹ��ͷ��칹������һ��
                    end
                    if mod(site_pos(1),20) == 0
                        site_col =  site_pos(1)/20;
                    else
                        site_col = fix(site_pos(1)/20)+1;
                    end
                    if mod(site_pos(2),20) == 0
                        site_row =  site_pos(2)/20;
                    else
                        site_row = fix(site_pos(2)/20)+1;
                    end
                    col_f = abs(site_col - q_col);
                    row_f = abs(site_row - q_row);
                    q_f = col_f + row_f;% ��ֱǽ��ˮƽǽ������Ϊ�ܵ�ǽ
                    if ~isFemto_het
                        q_f_matrix(:,:,b_) = reshape(q_f,roi_height_pixels,roi_width_pixels);
                    else
                        if strcmp(eNodeB_sites(b_).site_type, 'indoor')
                            % ֻ��indoor��վ�������ֵ���������͵Ļ�վΪ0
                            q_f_matrix(InH_row,InH_col,b_) = q_f;
                        end
                    end
                end
                
                for s_ = 1:length(eNodeB_sites(b_).sectors)
                    
                    % �������ص㵽��վ�ľ���С�����ھ�������,��ʱ���ھ���������ص㵽��վ�����һ��
                    d_2d_tmp = distance_matrix(:,:,b_);
                    d_2d_in(find(d_2d_tmp(:,:)<d_2d_in(:,:))) = d_2d_tmp(find(d_2d_tmp(:,:)<d_2d_in(:,:)))/2;
                    
                    %����pathloss����ģ�飬�õ�·�������LOS�ֲ�
                    [cell_pathloss_data(:,:,s_idx),LOS_compare(:,:,s_idx)] = eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model.pathloss(distance_matrix(:,:,b_),h_ut,d_2d_in);
                    
                    % ����͸��ı�ʾΪ��ά�������ʽ������pathlossһ������
                    if ~SYS_config.isDouble
                        penetration_loss(:,:,s_idx) = penetration;
                    else
                        num_first_sectors = networkMacroscopicPathlossMap.num_first_sectors;
                        if s_idx <= num_first_sectors
                            penetration_loss(:,:,s_idx) = penetration;
                        else
                            penetration_loss(:,:,s_idx) = penetration2;
                        end
                    end
                    if SYS_config.isFemto
                        % femto��ÿ��һ��ǽ��͸���5dB
                        % �Է��칹��˵���ұߵ�penetration_lossΪ0
                        % ���칹��˵���ұߵ�penetration_loss(:,:,�ڶ�ϵͳ��վ)��ʾ�����վ�����ڵĴ�ǽ��ģ�q_f_matrix(:,:,�ڶ�ϵͳ��վ)=0
                        penetration_loss(:,:,s_idx) = penetration_loss(:,:,s_idx) + q_f_matrix(:,:,b_) * 5;
                    end
                    
                    %�����վ�߶Ȳ����������������
                    num_first_sites = networkMacroscopicPathlossMap.num_first_sites;
                    switch SYS_config.scene_type
                        case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
                            if b_ <= num_first_sites
                                eNodeB_sites(b_).sectors(s_).tx_height = SYS_config.site_height;
                                eNodeB_sites(b_).sectors(s_).antenna.max_antenna_gain=SYS_config.antenna.max_antenna_gain;
                            else
                                eNodeB_sites(b_).sectors(s_).tx_height = SYS_config.site_height2;
                                eNodeB_sites(b_).sectors(s_).antenna.max_antenna_gain=SYS_config.antenna.max_antenna_gain2;
                            end
                        otherwise
                            eNodeB_sites(b_).sectors(s_).tx_height = SYS_config.site_height;
                            eNodeB_sites(b_).sectors(s_).antenna.max_antenna_gain=SYS_config.antenna.max_antenna_gain;
                    end
                    % ����ˮƽƫ��
                    switch eNodeB_sites(b_).site_type
                        case 'indoor'
                            angle_grid = (180/pi)*(atan2((position_grid_meters(:,2)-eNodeB_sites(b_).pos(2)),eNodeB_sites(b_).sectors(s_).tx_height-SYS_config.UE_height)); %��Ч����ϵ��ת90�ȣ����ҷ���x����ʾtheta��ǰ����y����ʾphi������phi�Դ�ֱ���·���Ϊ�ο�ϵ
                            angle_grid_ue=180-(180/pi)*(atan2((position_grid_meters(:,2)-eNodeB_sites(b_).pos(2)),position_grid_meters(:,1)-eNodeB_sites(b_).pos(1)));
                        otherwise
                            angle_grid = (180/pi)*(atan2((position_grid_meters(:,2)-eNodeB_sites(b_).pos(2)),(position_grid_meters(:,1)-eNodeB_sites(b_).pos(1))))-eNodeB_sites(b_).sectors(s_).azimuth;
                            angle_grid_ue=180-angle_grid;%UE��BS��ˮƽ�ǻ���
                    end
                    horizontal_angle_grid_ue   = reshape(angle_grid_ue,roi_height_pixels,roi_width_pixels);
                    horizontal_angle_grid_s_ue = utils.miscUtils.wrapTo359(horizontal_angle_grid_ue);
                    %ת��Ϊ(-180,180]
                    horizontal_angle_grid_s_ue=utils.miscUtils.wrapToAll180(horizontal_angle_grid_s_ue);
                    if strcmp(SYS_config.antenna.antenna_gain_pattern,'NRAntennaBeamforming')
                        
                        % ˮƽƫת�ǵ�Ч��-180~180
                        horizontal_angle_grid   = reshape(angle_grid,roi_height_pixels,roi_width_pixels);
                        horizontal_angle_grid_s = utils.miscUtils.wrapTo359(horizontal_angle_grid);
                        horizontal_angle_grid_s = utils.miscUtils.wrapToAll180(horizontal_angle_grid_s);
                        
                        % ��ֱƫ�ǽǶȼ���
                        switch eNodeB_sites(b_).site_type
                            case 'indoor'
                                vertical_angle_grid_el_ue=(180/pi)*(atan2(distance_matrix(:,:,b_),eNodeB_sites(b_).sectors(s_).tx_height - SYS_config.UE_height));
                                vertical_angle_grid_el = (180/pi)*(atan2( sqrt((SYS_config.site_height - SYS_config.UE_height).^2+(position_grid_meters(:,2)-eNodeB_sites(b_).pos(2)).^2),position_grid_meters(:,1)-eNodeB_sites(b_).pos(1)));%theta��ˮƽ����Ϊ�ο�ϵ
                                vertical_angle_grid_el   = reshape(vertical_angle_grid_el,roi_height_pixels,roi_width_pixels);
                            otherwise
                                vertical_angle_grid_el = (180/pi)*(atan2( distance_matrix(:,:,b_),h_ut(:,:)-eNodeB_sites(b_).sectors(s_).tx_height));
                                vertical_angle_grid_el_ue=(180/pi)*(atan2( distance_matrix(:,:,b_),eNodeB_sites(b_).sectors(s_).tx_height-h_ut(:,:)));
                        end
                        
                        phi_1=horizontal_angle_grid_s_ue+rnd_phi_panel1;
                        phi_2=horizontal_angle_grid_s_ue+rnd_phi_panel2;
                        %��������panel�����ת�Ǻ���Ҫ��֤���巽λ����-180��180��
                        phi_1 = utils.miscUtils.wrapToAll180(phi_1);
                        phi_2 = utils.miscUtils.wrapToAll180(phi_2);
                        UE_test = network_elements.UE;%���ڵ���element gain����
                        antennas.antenna.attach_antenna_to_UE(UE_test,SYS_config);
                        
                        switch eNodeB_sites(b_).site_type
                            case 'macro'
                                %��������panel��Ӧ��ˮƽ��ת�ǣ���ѡ��������������UEƫת��

                                sector_antenna_gain(:,:,s_idx) = eNodeB_sites(b_).sectors(s_).antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s, SYS_config.tilt);
                                UE_antenna_gain_1(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_1);
                                UE_antenna_gain_2(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_2);
                                UE_antenna_gain(:,:,s_idx)  = max(UE_antenna_gain_1(:,:,s_idx),UE_antenna_gain_2(:,:,s_idx));
                                
                                switch SYS_config.attatch_mode%���ý���ģʽ
                                    
                                    case 0 % ֻͨ��·�����
                                        sector_antenna_gain(:,:,s_idx)=0;
                                        UE_antenna_gain(:,:,s_idx)=0;
                                    case 1 %ͨ��·����BS�˵����ߵ�Ԫ�������
                                        UE_antenna_gain(:,:,s_idx)=0;
                                    case 2 %ͨ��·����UE�˵����ߵ�Ԫ�������
                                        sector_antenna_gain(:,:,s_idx)=0;
                                    case 3 %ͨ��·����BS��UE�����ߵ�Ԫ�������
                                end
                            case 'micro'
                                if ~SYS_config.isManhattan
                                    sector_antenna_gain(:,:,s_idx) = eNodeB_sites(b_).sectors(s_).antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s);
                                else
                                    sector_antenna_gain(:,:,s_idx) = 6;
                                end
                                UE_antenna_gain_1(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_1);
                                UE_antenna_gain_2(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_2);
                                UE_antenna_gain  = max(UE_antenna_gain_1,UE_antenna_gain_2);
                                
                                switch SYS_config.attatch_mode
                                    case 0
                                        sector_antenna_gain(:,:,s_idx)=0;
                                        UE_antenna_gain(:,:,s_idx)=0;
                                    case 1
                                        UE_antenna_gain(:,:,s_idx)=0;
                                    case 2
                                        sector_antenna_gain(:,:,s_idx)=0;
                                    case 3
                                end
                                
                            case 'indoor'
                                sector_antenna_gain(:,:,s_idx) = eNodeB_sites(b_).sectors(s_).antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s);
                                UE_antenna_gain_1(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_1);
                                UE_antenna_gain_2(:,:,s_idx)=UE_test.antenna.elementPattern(vertical_angle_grid_el_ue,phi_2);
                                
                                UE_antenna_gain  = max(UE_antenna_gain_1,UE_antenna_gain_2);
                                switch SYS_config.attatch_mode
                                    case 0
                                        sector_antenna_gain = zeros(size(sector_antenna_gain));
                                        UE_antenna_gain=zeros(size(UE_antenna_gain));
                                    case 1
                                        
                                        UE_antenna_gain=zeros(size(UE_antenna_gain));
                                    case 2
                                        sector_antenna_gain=zeros(size(sector_antenna_gain));
                                    case 3
                                end
                        end
                    else %���ʹ����������ģ��
                        error('"%s" antenna is not supported',SYS_config.antenna.antenna_gain_pattern);
                    end
                    
                    networkMacroscopicPathlossMap.sector_idx_mapping(eNodeB_sites(b_).sectors(s_).eNodeB_id,:) = [b_ s_];
                    networkMacroscopicPathlossMap.site_sector_mapping(b_,s_)  = eNodeB_sites(b_).sectors(s_).eNodeB_id;
                    
                    s_idx = s_idx + 1;
                end
            end
            
            %% ΢С�����ĵ������ص����
            micro_centre_distance_matrix = zeros(size(distance_matrix));
            for b_ = 1:length(eNodeB_sites)
                if strcmp(eNodeB_sites(b_).site_type,'micro')
                    if ~SYS_config.isManhattan
                        micro_centre__distances = sqrt((position_grid_meters(:,1)-eNodeB_sites(b_).parent_centre_pos(1)).^2 + (position_grid_meters(:,2)-eNodeB_sites(b_).parent_centre_pos(2)).^2);
                        micro_centre_distance_matrix(:,:,b_) = reshape(micro_centre__distances,roi_height_pixels,roi_width_pixels);% ΢С�����ĵ���ͼ�ϸ������ص�ľ���
                    end
                end
            end

            %% ��������
            cell_pathloss_data(isnan(cell_pathloss_data) | (cell_pathloss_data<0)) = 0;
            
            networkMacroscopicPathlossMap.micro_centre_distance_matrix = micro_centre_distance_matrix;% ΢С�����ĵ���ͼ�ϸ������ص�ľ���
            networkMacroscopicPathlossMap.distance_matrix = distance_matrix;% ��վ�������ص����
            networkMacroscopicPathlossMap.H_height = h_ut;% ���ص�߶�
            networkMacroscopicPathlossMap.H_orientation = rnd_phi_panel1;% ��ת��
            networkMacroscopicPathlossMap.penetration_loss = penetration_loss;% ���𣬺�·����һ����ģ�ľ���
            networkMacroscopicPathlossMap. pathloss_900 = cell_pathloss_data;% ���������������·��
            % �������·��
            networkMacroscopicPathlossMap.pathloss(:,:,eNodeB_id_set) = cell_pathloss_data - sector_antenna_gain-UE_antenna_gain+SYS_config.cable_loss;% �������������·��
            networkMacroscopicPathlossMap.LOS_compare(:,:,eNodeB_id_set) = LOS_compare;% LOS�ĸ��ʾ���
            networkMacroscopicPathlossMap.BS_antenna_gain = sector_antenna_gain;% ��վ����������
            networkMacroscopicPathlossMap.UE_antenna_gain = UE_antenna_gain;
        end
        
        function macroscopic_pathloss_model = generateMacroscopicPathlossModel(SYS_config,macroscopic_pathloss_model_name,frequency,macroscopic_pathloss_model_settings)
            % ͨ��macroscopic_pathloss_model_name�������ú��ʵ�·��ģ��
            print_output = true;
            switch macroscopic_pathloss_model_name
                case 'TS36942'% LTE·��ģ��
                    macroscopic_pathloss_model = macroscopic_pathloss_models.TS36942PathlossModel(frequency,macroscopic_pathloss_model_settings.environment);
                    if print_output && SYS_config.debug_level>=1
                        fprintf('TS 36.942-recommended pathloss model, %s environment\n',macroscopic_pathloss_model_settings.environment);
                    end
                case 'TS38900'% 5G·��ģ��
                    macroscopic_pathloss_model = macroscopic_pathloss_models.TS38900PathlossModel(frequency,macroscopic_pathloss_model_settings.environment,SYS_config);
                    if print_output && SYS_config.debug_level>=1
                        fprintf('TS 38.900-recommended pathloss model, %s environment\n',macroscopic_pathloss_model_settings.environment);
                    end
                case 'TS38901'% 5G·��ģ��
                    macroscopic_pathloss_model = macroscopic_pathloss_models.TS38901PathlossModel(frequency,macroscopic_pathloss_model_settings.environment,SYS_config);
                    if print_output && SYS_config.debug_level>=1
                        fprintf('TS 38.901-recommended pathloss model, %s environment\n',macroscopic_pathloss_model_settings.environment);
                    end 
                otherwise
                    error('"%s" macroscopic pathloss model not supported',macroscopic_pathloss_model_name);
            end
        end
    end
end
