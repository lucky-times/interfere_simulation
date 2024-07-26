classdef UE < handle
    % UEʵ���ļ�
    % ���幦�ܣ����ڼ���SINR

    properties

        id                    % UE��id
        pos                   % UE��λ��(x,y)
        attached_site         % UE������eNodeվַ
        attached_sector_idx   % UE������eNodeվַ��
        attached_eNodeB       % UE������eNodeB
        height                % UE�ĸ߶�
        orientation           % UE����ת��
        antenna               % UE������eNodeB
        walking_model         % UE���˶�ģ��
        downlink_channel      % UE�������ŵ�
        RB_grid               % UE��RB�����
        ul_ad_TX_power   %UE���й��غ�ķ��书��

        uplink_channel        % UE�������ŵ�
        PCe                   % ���й������
        receiver_noise_figure % UE�Ľ�������ϵ��
        BS_receiver_noise_figure % ����ʱBS���ն˵�����
        thermal_noise_W_RB    % ����������λW��
        BS_thermal_noise_W_RB
        UE_tx_power           % UE����ʱ�ķ��书��
        penetration_loss      % ��UE��������վ�Ĵ�͸���
        rx_antenna_gain       % UE������������
        beam_antenna_gain     % UE������������

        clock                 % ����ʱ��

        wideband_SINR        % �������м����SINR
        wideband_SINR2       % ��������������᷽�������SINR
        ul_wideband_SINR     % �������м����SINR
        ul_dl_wideband_SINR  % �������и�������ʱ��SINR
        dl_ul_wideband_SINR  % �������и�������ʱ��SINR
        bs_bs_pathloss       % ��վ�Ի�վ��������
        ue_ue_pathloss       % ue��ue��������

        dl_ul_UE_gain        %���и������еı�������·��UE���书��
        dl_ul_BS_gain        %���и������еı�������·��BS���书��

        ul_dl_UE_gain        %���и������еı�������·��UE���书��
        ul_dl_BS_gain        %���и������еı�������·��BS���书��


        SINR_TTI             % ����ÿ��TTI�е�����SINR
        ul_SINR_TTI          % ����ÿ��TTI�е�����SINR
        ul_dl_SINR_TTI       % ����ÿ��TTI�е����и�������SINR
        dl_ul_SINR_TTI       % ����ÿ��TTI�е����и�������SINR
        attached_eNodeB_id_TTI % ����ÿ��TTI��UE���ӵĻ�վ
        deactivate_UE        % ͳ�Ƹ�UE���
        cell_change          % �ṹ����С���л��й�

        interfering_UE
        interfering_BS

    end

    methods
        % ��ʼ�����캯��
        function obj = UE
            obj.deactivate_UE             = false;
            obj.cell_change.requested     = false;
            obj.cell_change.target_eNodeB = [];
        end

        function print(obj)
            if isempty(obj.attached_site)
                fprintf('User %d, (%d,%d), not attached to an eNodeB\n',obj.id,obj.pos(1),obj.pos(2));
            else
                fprintf('User %d, (%d,%d), Site %d, sector %d (eNodeB %d)\n',obj.id,obj.pos(1),obj.pos(2),obj.attached_site.id,obj.attached_sector_idx,obj.attached_eNodeB);
            end
            obj.walking_model.print;
        end

        %% ��ձ���
        function clear(obj)
            obj.attached_site             = [];
            obj.attached_eNodeB           = [];
            obj.walking_model             = [];
            obj.downlink_channel          = [];
            obj.RB_grid                   = [];
            obj.uplink_channel            = [];
            obj.clock                     = [];
        end

        %% �����ƶ�ģ�ͽ����ƶ�
        function move(obj)
            new_pos = obj.walking_model.move(obj.pos);
            obj.pos = new_pos;
        end

        function move_back(obj)
            old_pos = obj.walking_model.move_back(obj.pos);
            obj.pos = old_pos;
        end
        %% �ж�UE�Ƿ���ROI��
        function UE_in_roi = is_in_roi(obj,config,roi_x_range,roi_y_range)

            roi_x = obj.downlink_channel.macroscopic_pathloss_model.roi_x;
            roi_y = obj.downlink_channel.macroscopic_pathloss_model.roi_y;
            data_res = obj.downlink_channel.macroscopic_pathloss_model.data_res;
            sector_assignment = obj.downlink_channel.macroscopic_pathloss_model.sector_assignment;
            sector_assignment_double = obj.downlink_channel.macroscopic_pathloss_model.sector_assignment_double;
            num_first_UEs = obj.downlink_channel.macroscopic_pathloss_model.num_first_UEs;

            UE_pos_pix = NR_common_pos_to_pixel( obj.pos, [roi_x(1) roi_y(1)], data_res);

            try % ��������Խ��Ĵ����ж�UE��������ص��Ƿ���ROI��
                if config.isDouble
                    if UE_id<=num_first_UEs
                        s_idx = sector_assignment(UE_pos_pix(2),UE_pos_pix(1));
                    else
                        s_idx = sector_assignment_double(UE_pos_pix(2),UE_pos_pix(1));
                    end
                else
                    s_idx = sector_assignment(UE_pos_pix(2),UE_pos_pix(1));
                end
            catch
                UE_in_roi = false;
                return;
            end

            % ͨ�����ص㲻����ȫ�ж�UE����λ���Ƿ���ROI��,��inh��(120.4750,27.1949)���ж���ROI��
            UE_pos_x = obj.pos(1);
            UE_pos_y = obj.pos(2);
            if UE_pos_x<roi_x_range(1) || UE_pos_x>roi_x_range(2)
                UE_in_roi = false;
                return;
            end
            if UE_pos_y<roi_y_range(1) || UE_pos_y>roi_y_range(2)
                UE_in_roi = false;
                return;
            end

            if s_idx == -1 % UE��ROI�ϵ����ڻ�վ�ķ���Χ�ڣ���չʾ�õ�GUI�ϱ���Ϊ��ɫ�Ĳ���
                UE_in_roi = false;
                return;
            else
                UE_in_roi = true;
            end

            if UE_in_roi
                num_first_UEs = obj.downlink_channel.macroscopic_pathloss_model.num_first_UEs;
                [~,~,new_eNodeB_id] = obj.downlink_channel.macroscopic_pathloss_model.cell_assignment(config,obj.id,obj.pos,num_first_UEs);
                if obj.attached_eNodeB.eNodeB_id ~= new_eNodeB_id % ����������ȴ���Ҫ�����л�
                    obj.cell_change.requested = true;
                    obj.cell_change.target_eNodeB = new_eNodeB_id;
                end
            end
        end
        %% �л�����
        function start_handover(obj,new_eNodeB)
            obj.attached_eNodeB.deattachUser(obj);
            new_eNodeB.attachUser(obj);
        end

        %% ����ÿ��TTI��UE��SINR
        function save_SINR_in_TTI(obj,config)

            if config.asynchronization && config.isDouble
                obj.ul_dl_SINR_TTI(obj.clock.current_TTI) = obj.ul_dl_wideband_SINR;
                obj.dl_ul_SINR_TTI(obj.clock.current_TTI) = obj.dl_ul_wideband_SINR;
            else
                obj.SINR_TTI(obj.clock.current_TTI) = obj.wideband_SINR;
                obj.ul_SINR_TTI(obj.clock.current_TTI) = obj.ul_wideband_SINR;
            end

            obj.attached_eNodeB_id_TTI(obj.clock.current_TTI) = obj.attached_eNodeB.eNodeB_id;
        end

        %% ����4������������SINR
        % ����˼·��������񾶵�CL������ͬƵ���ž���CL��Ȼ������Ƶ���ž���CL��
        % ���ٸ��ݷ��书������շ����źŹ��ʣ������źŹ��ʣ��̶����SINR

        % ��������SINR
        function [IoT interference_level BS_Rx_power BS_Rx_intf_power ul_wideband_loop_SINR] = up_link_quality_model(obj,config,UEs,bs_beam_gain,ue_beam_gain)
            % ���������
            % config����������ϵͳ
            % UEs��UE�б�
            % bs_beam_gain��BS������������
            % ue_beam_gain��UE������������
            %
            % ���������
            % IoT��IoTֵ
            % interference_level�������źž���
            % ul_wideband_loop_SINR��UL����SINR
            interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors;% ���л�վ
            there_are_interferers            = ~isempty(interfering_eNodeBs);

            user_penetration_loss = obj.downlink_channel.macroscopic_penetration_loss; % �õ������źŵĴ���
            user_macroscopic_pathloss = obj.downlink_channel.macroscopic_pathloss + user_penetration_loss; % �õ������źŵ�·�𣨴�+��Ԫ����+����
            user_shadow_fading_loss = obj.downlink_channel.shadow_fading_pathloss; % �õ������źŵ���Ӱ˥��
            %�õ�UE��CL
            user_CL = user_macroscopic_pathloss + user_shadow_fading_loss-bs_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id)-ue_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id);
            user_CL_linear = 10^(0.1*user_CL);%ת��Ϊ����ֵ

            the_RB_grid  = obj.downlink_channel.RB_grid;
            nRB          = floor(the_RB_grid.n_RB/config.n_UE_served_per_BS);
            nSC          = nRB*2;

            % ��������
            thermal_noise_watts_per_half_RB = obj.BS_thermal_noise_W_RB/2;
            thermal_noise_watts = thermal_noise_watts_per_half_RB*nSC;

            % ÿ��UE������书��
            TX_power = config.UE_tx_power;
            %����ͬƵ���ż���Ƶ����
            if there_are_interferers
                % ����Ƶ����UE��ͬһ�������վ�µ�UE���ţ�
                interfering_ue_id = [];
                tmp2 = [obj.attached_eNodeB.attached_UEs_vector.id];
                index_in_cell = find(tmp2 == obj.id);% index_in_cellΪUE��������С����UE��˳��

                interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors;
                for s_ = 1:length(interfering_eNodeBs)
                    tmp = [interfering_eNodeBs(s_).attached_UEs_vector.id];
                    if length(tmp) == length(tmp2)% ��С�������С��UE�����
                        if isempty(interfering_ue_id)
                            interfering_ue_id = tmp(index_in_cell);
                        else
                            interfering_ue_id = [interfering_ue_id tmp(index_in_cell)];
                        end
                    else% ��С�������С��UE�������,ʵ���ϲ�֧��ÿ��С��UE������ͬ�����Σ�����ֻ�������һ�ֽ������

                        if isempty(interfering_ue_id)
                            interfering_ue_id = tmp(NR_randperm(length(tmp),1));
                        else
                            interfering_ue_id = [interfering_ue_id tmp(NR_randperm(length(tmp),1))];
                        end
                    end
                end
                interfering_ue_id = sort(interfering_ue_id); % Ϊ����UE����

                UE_pos_vector = obj.downlink_channel.macroscopic_pathloss_model.UE_pos_vector;
                interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);% ����UE��λ��

                % ����UE����ǰUE�ķ���BS��������
                interfering_macroscopic_pathloss_UE = obj.downlink_channel.ul_interfering_macroscopic_pathloss(obj.attached_eNodeB.eNodeB_id,interfering_ue_id);
                interfering_penetration_loss = obj.downlink_channel.ul_interfering_penetration_loss(obj.attached_eNodeB.eNodeB_id,interfering_ue_id);
                interfering_shadow_fading_loss = obj.downlink_channel.ul_interfering_shadow_fading_pathloss(obj.attached_eNodeB.parent_eNodeB.id,interfering_ue_id);
                interfering_CL = interfering_macroscopic_pathloss_UE + interfering_penetration_loss + interfering_shadow_fading_loss-bs_beam_gain(interfering_ue_id,obj.attached_eNodeB.eNodeB_id)'-ue_beam_gain(interfering_ue_id,obj.attached_eNodeB.eNodeB_id)';
                interfering_CL_linear = 10.^(0.1*interfering_CL);

                % ����UE�����Է���BS��������
                interfering_sector_id = fix((interfering_ue_id-1)/config.UE_per_eNodeB)+1;
                interfering_enb_id = obj.downlink_channel.macroscopic_pathloss_model.sector_idx_mapping(interfering_sector_id,1);
                interfering_enb_id = interfering_enb_id';
                for u_ = 1:length(interfering_ue_id)
                    interfering_pathloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_pathloss_eNodeB(interfering_ue_pos(u_,:),interfering_sector_id(u_));
                    interfering_penetrationloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_penetration_loss_eNodeB(interfering_ue_pos(u_,:),interfering_sector_id(u_));
                    interfering_shadow_fading_UE_to_sector(u_) = obj.downlink_channel.shadow_fading_model.get_pathloss(interfering_ue_pos(u_,:),interfering_enb_id(u_));
                    rx_antenna_gain(u_)=ue_beam_gain(interfering_ue_id(u_),UEs(interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);
                    beam_antenna_gain(u_)=bs_beam_gain(interfering_ue_id(u_),UEs(interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);
                end
                interfering_CL_UE_to_sector = interfering_pathloss_UE_to_sector + interfering_penetrationloss_UE_to_sector + interfering_shadow_fading_UE_to_sector-rx_antenna_gain-beam_antenna_gain;
                interfering_CL_UE_to_sector_linear = 10.^(0.1*interfering_CL_UE_to_sector);

                % ���ʿ���
                P_min = 10^(-4)/1000;
                R_min = P_min/TX_power;
                CL_x_ile = 88 + 10*log10(200/(config.bandwidth/1e6));
                CL_x_ile_linear = 10.^(0.1*CL_x_ile);
                Gama = config.Gama;
                TX_power_PC = TX_power * min(1,max(R_min,(user_CL_linear/CL_x_ile_linear)^Gama));
                interfering_TX_power_PC = TX_power .* min(1,max(R_min,(interfering_CL_UE_to_sector_linear./CL_x_ile_linear).^Gama));
                % ���빦�����
                TX_power_PC = TX_power_PC * 10.^(0.1*obj.PCe);
                for u_=1:length(interfering_TX_power_PC)
                    interfering_TX_power_PC(u_) = interfering_TX_power_PC(u_) * 10.^(0.1*UEs(interfering_ue_id(u_)).PCe);
                end

                obj.UE_tx_power = TX_power_PC;
                %������չ���
                RX_power = TX_power_PC/user_CL_linear;
                %ͬƵ���Ž��չ���
                interfering_RX_power = interfering_TX_power_PC./interfering_CL_linear;
                %����Ƶ���ţ�˫ϵͳ�£���һ��Ƶ�ʻ�վ�ĸ��ţ�
                if config.isDouble
                    ACIR_dB = config.ACIR_dB;
                    ACIR_linear = 10^(0.1*ACIR_dB);

                    % ����UE
                    % ΪӦ��һ����վͬʱ������UE�����⣬���ｫһ��С����ͬʱ�����UE��Ϊһ��
                    % ��index_in_cell=6��n_UE_served_per_BS=3����4,5,6Ϊһ��
                    if mod(index_in_cell,config.n_UE_served_per_BS) == 0
                        index_in_cell2 = index_in_cell-(config.n_UE_served_per_BS-1):index_in_cell;
                    else
                        yu = mod(index_in_cell,config.n_UE_served_per_BS);% yu=�����������UE�ڸ����ڵ�˳��
                        index_in_cell2 = index_in_cell-(yu-1):index_in_cell+(config.n_UE_served_per_BS-yu);
                    end

                    another_interfering_ue_id = [];
                    another_interfering_eNodeBs = obj.attached_eNodeB.ad_interf_eNodeB_sectors;
                    for s_ = 1:length(another_interfering_eNodeBs)
                        tmp = [another_interfering_eNodeBs(s_).attached_UEs_vector.id];
                        if length(tmp) == length(tmp2)% ��С�������С��UE�����
                            if isempty(another_interfering_ue_id)
                                another_interfering_ue_id = tmp(index_in_cell2);
                            else
                                another_interfering_ue_id = [another_interfering_ue_id tmp(index_in_cell2)];
                            end
                        else% ��С�������С��UE������ȣ�ʵ���ϲ�֧��ÿ��С��UE������ͬ�����Σ�����ֻ�������һ�ֽ������
                            if isempty(another_interfering_ue_id)
                                another_interfering_ue_id = tmp(NR_randperm(length(tmp),config.n_UE_served_per_BS));
                            else
                                another_interfering_ue_id = [another_interfering_ue_id tmp(NR_randperm(length(tmp),config.n_UE_served_per_BS))];
                            end
                        end
                    end
                    another_interfering_ue_id = sort(another_interfering_ue_id);
                    another_interfering_ue_pos = UE_pos_vector(another_interfering_ue_id,:);% ����UE��λ��

                    % ����UE����ǰUE�ķ���BS��������
                    another_interfering_macroscopic_pathloss_UE = obj.downlink_channel.ul_interfering_macroscopic_pathloss(obj.attached_eNodeB.eNodeB_id,another_interfering_ue_id);
                    another_interfering_penetration_loss = obj.downlink_channel.ul_interfering_penetration_loss(obj.attached_eNodeB.eNodeB_id,another_interfering_ue_id);
                    another_interfering_shadow_fading_loss = obj.downlink_channel.ul_interfering_shadow_fading_pathloss(obj.attached_eNodeB.parent_eNodeB.id,another_interfering_ue_id);
                    another_interfering_CL = another_interfering_macroscopic_pathloss_UE + another_interfering_penetration_loss + another_interfering_shadow_fading_loss-bs_beam_gain(another_interfering_ue_id,obj.attached_eNodeB.eNodeB_id)'-ue_beam_gain(another_interfering_ue_id,obj.attached_eNodeB.eNodeB_id)';
                    another_interfering_CL_linear = 10.^(0.1*another_interfering_CL);

                    % ����UE�����Է���BS��������
                    another_interfering_sector_id = fix((another_interfering_ue_id-1)/config.UE_per_eNodeB)+1;
                    another_interfering_enb_id = obj.downlink_channel.macroscopic_pathloss_model.sector_idx_mapping(another_interfering_sector_id,1);
                    another_interfering_enb_id = another_interfering_enb_id';

                    for u_ = 1:length(another_interfering_ue_id)
                        another_interfering_pathloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_pathloss_eNodeB(another_interfering_ue_pos(u_,:),another_interfering_sector_id(u_));
                        another_interfering_penetrationloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_penetration_loss_eNodeB(another_interfering_ue_pos(u_,:),another_interfering_sector_id(u_));
                        another_interfering_shadow_fading_UE_to_sector(u_) = obj.downlink_channel.shadow_fading_model.get_pathloss(another_interfering_ue_pos(u_,:),another_interfering_enb_id(u_));
                        rx_antenna_gain2(u_)=ue_beam_gain(another_interfering_ue_id(u_),UEs(another_interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);
                        beam_antenna_gain2(u_)=bs_beam_gain(another_interfering_ue_id(u_),UEs(another_interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);
                    end
                    another_interfering_CL_UE_to_sector = another_interfering_pathloss_UE_to_sector + another_interfering_penetrationloss_UE_to_sector + another_interfering_shadow_fading_UE_to_sector-rx_antenna_gain2-beam_antenna_gain2;
                    another_interfering_CL_UE_to_sector_linear = 10.^(0.1*another_interfering_CL_UE_to_sector);

                    % ���ʿ���
                    another_interfering_TX_power_PC = TX_power .* min(1,max(R_min,(another_interfering_CL_UE_to_sector_linear./CL_x_ile_linear).^Gama));
                    % ���빦�����
                    for u_=1:length(another_interfering_TX_power_PC)
                        another_interfering_TX_power_PC(u_) = another_interfering_TX_power_PC(u_) * 10.^(0.1*UEs(another_interfering_ue_id(u_)).PCe);
                    end
                    %��Ƶ���Ž��չ���
                    another_interfering_RX_power = another_interfering_TX_power_PC./another_interfering_CL_linear./ACIR_linear;

                    [~,n_i] = size(interfering_RX_power);
                    [~,n_ai] = size(another_interfering_RX_power);
                    if n_ai>n_i
                        % ����ϵͳ���Ż�վ���ȱ�����ϵͳ�Ķ�
                        interfering_RX_power(:,(n_i+1:n_ai)) = 0;
                    else
                        another_interfering_RX_power(:,(n_ai+1:n_i)) = 0;
                    end

                    switch config.interference_type
                        case 0
                            total_interfering_power = interfering_RX_power;
                        case 1
                            total_interfering_power = another_interfering_RX_power;
                        case 2
                            total_interfering_power = interfering_RX_power + another_interfering_RX_power;
                    end
                    % ����ACIR�Ƿ�ʹ�ò���ѭ�����е���
                    if config.isACIR_loop
                        loop_ACIR_dB = config.ACIR_lower:config.loop_step:config.ACIR_upper;
                        loop_ACIR_linear = 10.^(0.1*loop_ACIR_dB);
                        another_interfering_RX_power_1 = another_interfering_RX_power * ACIR_linear; %����ACIR����aggressorϵͳ���ĸ���ǿ��
                        for i=1:length(loop_ACIR_dB)
                            loop_another_interfering_RX_power(i,:) = another_interfering_RX_power_1 / loop_ACIR_linear(i);
                            loop_total_interfering_power(i,:) = interfering_RX_power + loop_another_interfering_RX_power(i,:);
                            ul_wideband_loop_SINR(i,:) = 10*log10(RX_power/(sum(loop_total_interfering_power(i,:))+thermal_noise_watts));
                        end
                        %ֻ��ͬƵ����ʱ��SINR
                        ul_wideband_loop_SINR(length(loop_ACIR_dB)+1,:) = 10*log10(RX_power/(sum(interfering_RX_power(:))+thermal_noise_watts));
                    else
                        ul_wideband_loop_SINR = 0;
                    end

                    % �����ź�
                    another_interfering_RX_power_2 = another_interfering_RX_power * ACIR_linear;
                    interference_level = sum(another_interfering_RX_power_2(:));
                else
                    total_interfering_power = interfering_RX_power;
                    ul_wideband_loop_SINR = 0;
                    interference_level = 0;
                end

                % SINR
                obj.ul_wideband_SINR = 10*log10(RX_power/(sum(total_interfering_power(:))+thermal_noise_watts));
                IoT = 10*log10(sum(total_interfering_power(:))/thermal_noise_watts);% ����ÿ��С����10��UE��ѭ����UE1ʱ����BS1���IOT��I����UE11��21��31�ĸ���
                % Ȼ��ѭ��UE1��UE10����BS1ȡ10�ο��յõ��ĸ���

                BS_Rx_power = RX_power;
                BS_Rx_intf_power = sum(total_interfering_power(:));
            else
                RX_power = TX_power/user_CL_linear;
                obj.ul_wideband_SINR = 10*log10(RX_power/thermal_noise_watts);
                ul_wideband_loop_SINR = 0;

                BS_Rx_power = RX_power;
                BS_Rx_intf_power = 0;
            end
        end

        function dummy_up_link_quality_model(obj)
            obj.ul_wideband_SINR = NaN;
            obj.UE_tx_power = NaN;
        end

        % ����SINR
        function [signal_CL interfering_CL another_interfering_CL UE_Rx_power UE_Rx_intf_power wideband_loop_SINR] = down_link_quality_model(obj,config,bs_beam_gain,ue_beam_gain, SYS_config)
            % ���������
            % config����������ϵͳ
            % UEs��UE�б�
            % bs_beam_gain��BS������������
            % ue_beam_gain��UE������������
            %
            % ���������
            % signal_CL�������ź�CL
            % interfering_CL��victimϵͳ�����ź�CL
            % another_interfering_CL��aggressorϵͳ�����ź�CL
            % wideband_loop_SINR��UE����SINR

            interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors;% �õ�victimϵͳ����С��
            another_interfering_eNodeBs = obj.attached_eNodeB.ad_interf_eNodeB_sectors; % �õ�aggressorϵͳ����С��
            there_are_interferers = ~isempty(interfering_eNodeBs);

            obj.penetration_loss = obj.downlink_channel.macroscopic_penetration_loss; % ����

            d2d = pdist2(obj.attached_site.pos, obj.pos,'Euclidean');
            d3d = sqrt(d2d^2 + (obj.height - obj.attached_eNodeB.tx_height)^2);
            distance = d3d;
            frequency = 4.9e9/1e9;
            pl = 28 + 22*log10(distance) + 20*log10(frequency);
            vertical_angle_grid_el = (180/pi)*(atan2( d2d,obj.height-obj.attached_eNodeB.tx_height));
            vertical_angle_grid_el_ue=(180/pi)*(atan2( d2d,obj.attached_eNodeB.tx_height-obj.height));
            alpha = 90 - vertical_angle_grid_el;
            d3d_cos = d2d/cosd(alpha);
            pl2 = 28 + 22*log10(d3d_cos) + 20*log10(frequency);

            angle_grid = (180/pi)*(atan2((obj.pos(2)-obj.attached_site.pos(2)),(obj.pos(1)-obj.attached_site.pos(1))))-obj.attached_eNodeB.azimuth;
            horizontal_angle_grid = utils.miscUtils.wrapTo359(angle_grid);
            horizontal_angle_grid_s = utils.miscUtils.wrapTo180(horizontal_angle_grid);
            
            angle_grid_ue=180-angle_grid;%UE��BS��ˮƽ�ǻ���
            horizontal_angle_grid_s_ue = utils.miscUtils.wrapTo359(angle_grid_ue);
            horizontal_angle_grid_s_ue=utils.miscUtils.wrapToAll180(horizontal_angle_grid_s_ue);
            rnd_phi_panel1=180-360*rand(size(obj.attached_eNodeB.tx_height)); %����-180~180�������ת��λ
            rnd_phi_panel2=rnd_phi_panel1-180;%����panel��180��
            phi_1=horizontal_angle_grid_s_ue+rnd_phi_panel1;
            phi_2=horizontal_angle_grid_s_ue+rnd_phi_panel2;
            %��������panel�����ת�Ǻ���Ҫ��֤���巽λ����-180��180��
            phi_1 = utils.miscUtils.wrapToAll180(phi_1);
            phi_2 = utils.miscUtils.wrapToAll180(phi_2);


            BS_element_pattern = obj.attached_eNodeB.antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s, SYS_config.tilt);
            UE_antenna_gain_1=obj.antenna.elementPattern(vertical_angle_grid_el_ue,phi_1);
            UE_antenna_gain_2=obj.antenna.elementPattern(vertical_angle_grid_el_ue,phi_2);
            UE_antenna_gain= max(UE_antenna_gain_1,UE_antenna_gain_2);
            macroscopic_pathloss = pl - BS_element_pattern - UE_antenna_gain;
            macroscopic_pathloss2 = pl2 - BS_element_pattern - UE_antenna_gain;
            user_macroscopic_pathloss        = macroscopic_pathloss ;   % �������
            user_macroscopic_pathloss_linear = 10^(0.1*user_macroscopic_pathloss);% ת��Ϊ����ֵ
            user_macroscopic_pathloss_linear2 = 10^(0.1*macroscopic_pathloss2);% ת��Ϊ����ֵ
            user_shadow_fading_loss          = obj.downlink_channel.shadow_fading_pathloss;% ��Ӱ˥��
            user_shadow_fading_loss = 0;
            user_shadow_fading_loss_linear   = 10^(0.1*user_shadow_fading_loss); % ת��Ϊ����ֵ

            signal_CL = user_macroscopic_pathloss + user_shadow_fading_loss; % �õ������źŵ�CL

            the_RB_grid  = obj.downlink_channel.RB_grid;
            nRB          = floor(the_RB_grid.n_RB/config.n_UE_served_per_BS);
            nSC          = nRB*2;

            % ��������
            thermal_noise_watts_per_half_RB = obj.thermal_noise_W_RB/2;
            thermal_noise_watts = thermal_noise_watts_per_half_RB*nSC;
            % ���书��
            TX_power = obj.attached_eNodeB.max_power/config.n_UE_served_per_BS;
            RX_power = TX_power / user_macroscopic_pathloss_linear / user_shadow_fading_loss_linear;
            RX_power2 = TX_power / user_macroscopic_pathloss_linear2;

            if there_are_interferers % ��������ź�CL
                parent_sites                            = [interfering_eNodeBs.parent_eNodeB];
                parent_sites_id                         = [parent_sites.id];
                interfering_eNodeB_ids                  = [interfering_eNodeBs.eNodeB_id];
                interfering_RB_grids                    = [interfering_eNodeBs.RB_grid];
                interfering_power_allocations_data      = [interfering_RB_grids.power_allocation];
                interfering_power_allocations_signaling = [interfering_RB_grids.power_allocation_signaling];
                interfering_power_allocations = interfering_power_allocations_data + interfering_power_allocations_signaling;% n_����enb�У�n_RB��
                interfering_power_allocations = sum(interfering_power_allocations);
                for i = 1:length(parent_sites)
                    interfer_sites_pos(i, :) = parent_sites(i).pos;
                    inteNodeB_azimuth(i) = interfering_eNodeBs(i).azimuth;
                end
                d2d_het = pdist2(interfer_sites_pos, obj.pos,'Euclidean');
                d3d_het = sqrt(d2d_het.^2 + (obj.height - obj.attached_eNodeB.tx_height).^2);
                distance_het = d3d_het;
                pl_het = 28 + 22*log10(distance_het) + 20*log10(frequency);
%                 pl_het(1:2) =[150,150];
                vertical_angle_grid_el = (180/pi)*(atan2( d2d_het,obj.height-obj.attached_eNodeB.tx_height));
                vertical_angle_grid_el_ue=(180/pi)*(atan2( d2d_het,obj.attached_eNodeB.tx_height-obj.height));
%                 alpha_het = 90 - vertical_angle_grid_el;
%                 d3d_cos_het = d2d_het./cosd(alpha_het);
                n_i = d2d_het / d2d;

                pl2_het = 28 + 22*log10(d3d_cos) + 22*log10(cosd(alpha)./cosd(alpha./n_i) .* n_i)+ 20*log10(frequency);
                if sum(pl2_het- pl_het)>10
                    pl2_het- pl_het;
                end
                angle_grid = (180/pi)*(atan2((obj.pos(2)-interfer_sites_pos(:,2)),(obj.pos(1)-interfer_sites_pos(:,1))))-inteNodeB_azimuth';
%                 phi = angle_grid;
%                 phi(1) = phi(1)+120;
%                 phi(2) = phi(2)+240;
%                 for i = 3:length(phi)
%                     phi(i) = phi(i) + mod(i,3) * 120;
%                 end
                horizontal_angle_grid = utils.miscUtils.wrapTo359(angle_grid);
                horizontal_angle_grid_s = utils.miscUtils.wrapTo180(horizontal_angle_grid);
                
                angle_grid_ue=180-angle_grid;%UE��BS��ˮƽ�ǻ���
                horizontal_angle_grid_s_ue = utils.miscUtils.wrapTo359(angle_grid_ue);
                horizontal_angle_grid_s_ue=utils.miscUtils.wrapToAll180(horizontal_angle_grid_s_ue);
                rnd_phi_panel1=180-360*rand(length(parent_sites), 1); %����-180~180�������ת��λ
                rnd_phi_panel2=rnd_phi_panel1-180;%����panel��180��
                phi_1=horizontal_angle_grid_s_ue+rnd_phi_panel1;
                phi_2=horizontal_angle_grid_s_ue+rnd_phi_panel2;
                %��������panel�����ת�Ǻ���Ҫ��֤���巽λ����-180��180��
                phi_1 = utils.miscUtils.wrapToAll180(phi_1);
                phi_2 = utils.miscUtils.wrapToAll180(phi_2);


                BS_element_pattern = obj.attached_eNodeB.antenna.elementPattern(vertical_angle_grid_el, horizontal_angle_grid_s, SYS_config.tilt);
                BS_element_pattern((90-vertical_angle_grid_el)>(SYS_config.tilt+24)) = -inf;
                BS_element_pattern((90-vertical_angle_grid_el)<SYS_config.tilt) = -inf;
                UE_antenna_gain_1=obj.antenna.elementPattern(vertical_angle_grid_el_ue,phi_1);
                UE_antenna_gain_2=obj.antenna.elementPattern(vertical_angle_grid_el_ue,phi_2);
                UE_antenna_gain= max(UE_antenna_gain_1,UE_antenna_gain_2);
                
                macroscopic_pathloss_het = pl_het - BS_element_pattern - UE_antenna_gain;
%                 macroscopic_pathloss_het = macroscopic_pathloss_het(3:3:end);
                macroscopic_pathloss2_het = pl2_het - BS_element_pattern - UE_antenna_gain;
%                 macroscopic_pathloss2_het = real(macroscopic_pathloss2_het);
%                 macroscopic_pathloss2_het = macroscopic_pathloss2_het(3:3:end);

                interfering_macroscopic_pathloss_eNodeB        = macroscopic_pathloss_het ;%����·����λ��ʱ��
                interfering_shadow_fading_loss                 = obj.downlink_channel.interfering_shadow_fading_pathloss(parent_sites_id);
                interfering_shadow_fading_loss = 0;
                interfering_CL = interfering_macroscopic_pathloss_eNodeB + interfering_shadow_fading_loss;

                interfering_macroscopic_pathloss_eNodeB_linear = 10.^(0.1*interfering_macroscopic_pathloss_eNodeB);
                interfering_macroscopic_pathloss_eNodeB_linear2 = 10.^(0.1*macroscopic_pathloss2_het);
                interfering_shadow_fading_loss_linear          = 10.^(0.1*interfering_shadow_fading_loss);
                interfering_power = interfering_power_allocations ./ interfering_macroscopic_pathloss_eNodeB_linear';
                interfering_power2 = interfering_power_allocations ./ interfering_macroscopic_pathloss_eNodeB_linear2';
                if config.isDouble %�����˫ϵͳ������aggressorϵͳ�ĸ����ź�
                    another_parent_sites                            = [another_interfering_eNodeBs.parent_eNodeB];
                    another_parent_sites_id                         = [another_parent_sites.id];
                    another_interfering_eNodeB_ids                  = [another_interfering_eNodeBs.eNodeB_id];
                    another_interfering_RB_grids                    = [another_interfering_eNodeBs.RB_grid];
                    another_interfering_power_allocations_data      = [another_interfering_RB_grids.power_allocation];
                    another_interfering_power_allocations_signaling = [another_interfering_RB_grids.power_allocation_signaling];
                    another_interfering_power_allocations = another_interfering_power_allocations_data + another_interfering_power_allocations_signaling;
                    another_interfering_power = sum(another_interfering_power_allocations);

                    another_interfering_macroscopic_pathloss_eNodeB        = obj.downlink_channel.interfering_macroscopic_pathloss(another_interfering_eNodeB_ids) + obj.downlink_channel.interfering_penetration_loss(another_interfering_eNodeB_ids) - bs_beam_gain(obj.id,another_interfering_eNodeB_ids)'-ue_beam_gain(obj.id,another_interfering_eNodeB_ids)';
                    another_interfering_shadow_fading_loss                 = obj.downlink_channel.interfering_shadow_fading_pathloss(another_parent_sites_id);

                    another_interfering_CL = another_interfering_macroscopic_pathloss_eNodeB + another_interfering_shadow_fading_loss;
                    % Ϊ��ͬƵ�Աȣ�ɾ����UE�ķ���BSͬһ��С����BS
                    [~,index_min] = sort(another_interfering_CL);
                    another_interfering_CL(index_min(1)) = [];

                    another_interfering_macroscopic_pathloss_eNodeB_linear = 10.^(0.1*another_interfering_macroscopic_pathloss_eNodeB);
                    another_interfering_shadow_fading_loss_linear          = 10.^(0.1*another_interfering_shadow_fading_loss);

                    another_interfering_power = another_interfering_power ./ another_interfering_macroscopic_pathloss_eNodeB_linear' ./ another_interfering_shadow_fading_loss_linear';

                    [~,n_i] = size(interfering_power);
                    [~,n_ai] = size(another_interfering_power);
                    if n_ai>n_i
                        % ����ϵͳ���Ż�վ���ȱ�����ϵͳ�Ķ�
                        interfering_power(:,(n_i+1:n_ai)) = 0;
                    else
                        another_interfering_power(:,(n_ai+1:n_i)) = 0;
                    end




                    if config.isACIR_loop
                        loop_ACIR_dB = config.ACIR_lower:config.loop_step:config.ACIR_upper;
                        loop_ACIR_linear = 10.^(0.1*loop_ACIR_dB);
                        for i=1:length(loop_ACIR_dB)
                            loop_another_interfering_power(i,:) = another_interfering_power / loop_ACIR_linear(i);
                            loop_total_interfering_power(i,:) = interfering_power + loop_another_interfering_power(i,:);
                            wideband_loop_SINR(i,:) = 10*log10(RX_power/(sum(loop_total_interfering_power(i,:))+thermal_noise_watts));
                        end
                        wideband_loop_SINR(length(loop_ACIR_dB)+1,:) = 10*log10(RX_power/(sum(interfering_power(:))+thermal_noise_watts));
                    else
                        wideband_loop_SINR = 0;
                    end

                    ACIR_dB = config.ACIR_dB;
                    ACIR_linear = 10^(0.1*ACIR_dB);
                    another_interfering_power = another_interfering_power / ACIR_linear;

                    switch config.interference_type
                        case 0
                            total_interfering_power = interfering_power;% ֻ�б�ϵͳ����
                            total_interfering_power2 = interfering_power2;% ֻ�б�ϵͳ����
                        case 1
                            total_interfering_power = another_interfering_power;% ֻ�еڶ�ϵͳ����
                        case 2
                            total_interfering_power = interfering_power + another_interfering_power; %����ϵͳ�ĸ���
                    end
                else
                    total_interfering_power = interfering_power;% ֻ�б�ϵͳ����
                    total_interfering_power2 = interfering_power2;% ֻ�б�ϵͳ����
                    another_interfering_CL = zeros(size(interfering_power));
                    wideband_loop_SINR = 0;
                end

                % SINR
                obj.wideband_SINR = 10*log10(RX_power/(sum(total_interfering_power(:))+thermal_noise_watts));
                obj.wideband_SINR2 = 10*log10(RX_power2/(sum(total_interfering_power2(:))+thermal_noise_watts));
                UE_Rx_power = RX_power;
                UE_Rx_intf_power = sum(total_interfering_power(:));
            else
                wideband_loop_SINR = 0;
                obj.wideband_SINR = 10*log10(RX_power/thermal_noise_watts);

                UE_Rx_power = RX_power;
                UE_Rx_intf_power = 0;
            end

        end

        function dummy_down_link_quality_model(obj)
            obj.wideband_SINR = NaN;
        end

        %% ���и������У�UE����UE��ֻ��˫ϵͳʱ�Ż�ʹ��
        function [user_macroscopic_pathloss,ul_down_wideband_loop_SINR] = up_down_model(obj,config,UEs,bs_beam_gain,ue_beam_gain,networkPathlossMap)
            interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors; %�õ�����eNodeB
            obj.penetration_loss = obj.downlink_channel.macroscopic_penetration_loss; %�õ���UE�������BS�Ĵ���
            user_macroscopic_pathloss        = obj.downlink_channel.macroscopic_pathloss_900;   % ��+��Ԫ����
            user_shadow_fading_loss          = obj.downlink_channel.shadow_fading_pathloss; %��Ӱ˥��
            user_CL = user_macroscopic_pathloss + user_shadow_fading_loss+ obj.penetration_loss -bs_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id)- ue_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id);%�����CL
            user_CL_linear = 10^(0.1*user_CL);%ת��Ϊ����ֵ

            the_RB_grid  = obj.downlink_channel.RB_grid;
            nRB          = floor(the_RB_grid.n_RB/config.n_UE_served_per_BS);
            nSC          = nRB*2;

            % ����
            thermal_noise_watts_per_half_RB = obj.thermal_noise_W_RB/2;
            thermal_noise_watts = thermal_noise_watts_per_half_RB*nSC;
            % �õ���վ���书��
            TX_power = obj.attached_eNodeB.max_power/config.n_UE_served_per_BS;
            % �õ���UE���յ�����BS�Ĺ���
            RX_power = TX_power /user_CL_linear;

            %��������Ż�վ�ķ��书��
            parent_sites                            = [interfering_eNodeBs.parent_eNodeB];
            parent_sites_id                         = [parent_sites.id];
            interfering_eNodeB_ids                  = [interfering_eNodeBs.eNodeB_id];
            interfering_RB_grids                    = [interfering_eNodeBs.RB_grid];
            interfering_power_allocations_data      = [interfering_RB_grids.power_allocation];
            interfering_power_allocations_signaling = [interfering_RB_grids.power_allocation_signaling];
            interfering_power_allocations = interfering_power_allocations_data + interfering_power_allocations_signaling;% n_����enb��n_RB��
            interfering_power = sum(interfering_power_allocations);

            interfering_macroscopic_pathloss_eNodeB        = obj.downlink_channel.interfering_macroscopic_pathloss(interfering_eNodeB_ids); %�������+��Ԫ����
            interfering_shadow_fading_loss                 = obj.downlink_channel.interfering_shadow_fading_pathloss(parent_sites_id);%��Ӱ˥��
            interfering_CL=interfering_macroscopic_pathloss_eNodeB+interfering_shadow_fading_loss+obj.downlink_channel.interfering_penetration_loss(interfering_eNodeB_ids) - bs_beam_gain(obj.id,interfering_eNodeB_ids)'-ue_beam_gain(obj.id,interfering_eNodeB_ids)';%������
            %interfering_PL =    interfering_macroscopic_pathloss_eNodeB-
            interfering_CL_linear = 10.^(0.1*interfering_CL);
            interfering_power = interfering_power ./interfering_CL_linear';%������յ��ĸ���BS����

            %% ��ϵͳ�ڲ������и��ŵ���ʽ������Ҫ���ݸ��ݸ���UE�뱻����UE��������ߵĵ�Ԫ���棬·����ģ�����ȫΪNLOS�������Ҽ�������Ӱ˥���Ӱ��
            another_interfering_eNodeBs = obj.attached_eNodeB.ad_interf_eNodeB_sectors; %�ҳ���ϵͳ�ĸ���BS
            another_interfering_eNodeB_ids = [another_interfering_eNodeBs.eNodeB_id]; %�ҵ��ڶ���ϵͳ�ĸ�������id��
            ser_group_num = mod(obj.id,config.UE_per_eNodeB); %�ҳ��������
            if ser_group_num==0
                ser_group_num=config.UE_per_eNodeB;
            end
            another_interfering_UE_ids = (another_interfering_eNodeB_ids-1).* config.UE_per_eNodeB + ser_group_num;%�ҳ�����UE��id��
            for u_=1:length(another_interfering_UE_ids)
                %�����UE������ʱ����ϵͳС��u_�ĸ�����
                intf_group =floor((another_interfering_UE_ids(u_)-1)/config.n_UE_served_per_BS)*config.n_UE_served_per_BS+1 :floor((another_interfering_UE_ids(u_)-1)/config.n_UE_served_per_BS)*config.n_UE_served_per_BS+config.n_UE_served_per_BS;
                %����ÿ������UE����ڸ�UE��ˮƽ�ǣ���ת��Ϊ(-180,180]
                for j=1:config.n_UE_served_per_BS     %j��ʾ��ϵͳС��������ĵ�j��UE
                    phi_up_down(u_,j)=atan2(UEs(intf_group(j)).pos(2)-obj.pos(2),UEs(intf_group(j)).pos(1)-obj.pos(1))./pi*180-obj.orientation;
                    phi_up_down(u_,j)=utils.miscUtils.wrapToAll180(phi_up_down(u_,j));
                    %����ÿ������UE����ڸ�UE�Ĵ�ֱ��
                    theta_up_down(u_,j) = atan2(sqrt((obj.pos(1)-UEs(intf_group(j)).pos(1)).^2+(obj.pos(2)-UEs(intf_group(j)).pos(2))^2),UEs(intf_group(j)).height-obj.height)./pi*180;
                    %����ÿ������UE����ڸ�UE�ĵ�Ԫ����
                    element_gain_interference(u_,j) = obj.antenna.elementPattern(theta_up_down(u_,j),phi_up_down(u_,j));
                    %�����UE�����ÿ������UE�Ĵ�ֱ����ˮƽ��
                    theta_opposite_observation(u_,j)=atan2(sqrt((obj.pos(1)-UEs(intf_group(j)).pos(1)).^2+(obj.pos(2)-UEs(intf_group(j)).pos(2))^2),obj.height-UEs(intf_group(j)).height)./pi*180;
                    phi_opposite_observation(u_,j)=atan2(obj.pos(2)-UEs(intf_group(j)).pos(2),obj.pos(1)-UEs(intf_group(j)).pos(1))./pi*180-UEs(intf_group(j)).orientation;
                    %% �Ƕȵ���Ϊ(-180,180]
                    phi_opposite_observation(u_,j)=utils.miscUtils.wrapToAll180(phi_opposite_observation(u_,j));

                    %�����UE�����ÿ������UE�ĵ�Ԫ����
                    element_gain_UE(u_,j) = obj.antenna.elementPattern(theta_opposite_observation(u_,j),phi_opposite_observation(u_,j));
                    %������ߵ�3d����
                    d_3d(u_,j) = sqrt((obj.pos(1)-UEs(intf_group(j)).pos(1)).^2+(obj.pos(2)-UEs(intf_group(j)).pos(2)).^2+(obj.height-UEs(intf_group(j)).height).^2);
                    %����pathloss
                    if obj.id <= networkPathlossMap.num_first_UEs
                        frequency=config.frequency2;
                    else
                        frequency=config.frequency;
                    end
                    switch config.scene_type %������NLOS·��ģ��
                        case {'InH','InH2','InH3'}%�����indoor�������������·��ģ��
                            shadow_fading_std = 8.03;%��������NLOSʱ����Ӱ˥��
                            %                             shadow_fading_std = 3.0;%��������LOSʱ����Ӱ˥��
                            shadow_fading_dB(u_,j)=normrnd(0,shadow_fading_std);%����������̬�ֲ�����Ӱ˥��
                            pl(u_,j)=17.3+38.3*log10(d_3d(u_,j))+24.9*log10(frequency/1e9)+shadow_fading_dB(u_,j);%NLOS
                            %                             pl(u_,j)=32.4+17.3*log10(d_3d(u_))+20*log10(frequency/1e9);%LOS
                        otherwise %XIAģ���в��� ���ߵ���ƽ������ˮƽ��ģ��,������Ӱ��������ĺ�������ġ���A Simplified Analytical Model for Predicting Path Loss in Urban and Suburban Environments��
                            h_roof = 12;%�ݶ��߶�12m
                            w=30;%�ֵ����30m
                            d=80;%���������ľ���
                            R=d_3d(u_,j); %���������ջ����
                            delta_hb= obj.height - h_roof;%���ջ�����ڷ����ĸ߶ȣ��ȷ�������Ϊ����
                            delta_hm = h_roof - UEs(intf_group(j)).height;%���������ڷ����ĸ߶Ȳ�
                            phi=-atan2(delta_hb,d);%����һ�Ž������Ե�������
                            lamda = 3e8/config.frequency; %����
                            x=w/2;%����UE���������ߵ�ˮƽ����
                            theta=atan2(delta_hm,x);
                            r=sqrt(delta_hm^2+x^2);%
                            pl(u_,j)=-10*log10((lamda/(2*sqrt(2)*pi*R))^2)-10*log10(lamda/(2*pi^2*r)*((1/theta)-1/(2*pi+theta))^2)-10*log((d/(2*pi*(R-d)))^2*lamda/sqrt(delta_hb^2+d^2)*(1/phi-1/(2*pi+phi))^2);
                    end

                    %���㴩͸���
                    switch config.scene_type
                        case {'UMA','UMI','RMa','UMa_to_UMi'}
                            penetration(u_,j)=obj.downlink_channel.macroscopic_penetration_loss+UEs(intf_group(j)).downlink_channel.macroscopic_penetration_loss;
                        case {'InH','InH2'}
                            penetration(u_,j)=0;
                        case {'UMi_to_InH','UMa_to_InH'} %�����������ͬһ���ڵ��û��޴�͸��ģ��������Urban�����µ��û�һ���ᴩһ��ǽ
                            L_glass = 2+0.2*config.frequency/1e9;
                            L_concrete = 5+4*config.frequency/1e9;
                            N_low = 0+4.4*randn(1);
                            low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_low;
                            L_IRRglass = 23+0.3*config.frequency/1e9;
                            N_high = 0+6.5*randn(1);
                            high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_high;
                            if rand(1)>0.5
                                penetration(u_,j)=low;
                            else
                                penetration(u_,j)=high;
                            end
                        case {'InH3'} %�����������ͬһ���ڵ��û��޴�͸��ģ��������Urban�����µ��û�һ���ᴩһ��ǽ
                            L_glass = 2+0.2*config.frequency/1e9;
                            L_concrete = 5+4*config.frequency/1e9;
                            N_low = 0+4.4*randn(1);
                            low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_low;
                            L_IRRglass = 23+0.3*config.frequency/1e9;
                            N_high = 0+6.5*randn(1);
                            high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_high;
                            if rand(1)>0.5
                                penetration(u_,j)=2*low;
                            else
                                penetration(u_,j)=2*high;
                            end

                    end
                    %����������

                    another_interfering_CL(u_,j) = pl(u_,j) + penetration(u_,j)-element_gain_interference(u_,j) - element_gain_UE(u_,j) - ue_beam_gain(obj.id,UEs(intf_group(j)).attached_eNodeB.eNodeB_id) - bs_beam_gain(obj.id,UEs(intf_group(j)).attached_eNodeB.eNodeB_id);
                    another_interfering_CL_linear(u_,j) = 10.^(0.1*another_interfering_CL(u_,j)); %ת��Ϊ����ֵ
                    %�������UE�Ĺ��غ��书�ʣ�ͨ��UE�Է����վ��CL��ȷ����
                    another_interfering_pathloss_UE_to_sector(u_,j) = obj.downlink_channel.macroscopic_pathloss_model.get_pathloss_eNodeB(UEs(intf_group(j)).pos,UEs(intf_group(j)).attached_eNodeB.eNodeB_id);%����UE�������С���Ĵ������
                    another_interfering_penetrationloss_UE_to_sector(u_,j) = obj.downlink_channel.macroscopic_pathloss_model.get_penetration_loss_eNodeB(UEs(intf_group(j)).pos,UEs(intf_group(j)).attached_eNodeB.eNodeB_id);%����UE�������С���Ĵ�͸���
                    another_interfering_shadow_fading_UE_to_sector(u_,j) = obj.downlink_channel.shadow_fading_model.get_pathloss(UEs(intf_group(j)).pos,UEs(intf_group(j)).attached_eNodeB.parent_eNodeB.id);%����UE�������С������Ӱ˥��
                    tx_antenna_gain(u_,j)=ue_beam_gain(intf_group(j),UEs(intf_group(j)).attached_eNodeB.eNodeB_id); %����UE�������С����UE�ನ����������
                    rx_antenna_gain(u_,j)=bs_beam_gain(intf_group(j),UEs(intf_group(j)).attached_eNodeB.eNodeB_id); %����UE�������С����BS�ನ����������


                    %����UE�������С���������ģ����ڼ��㹦�أ�
                    another_interfering_CL_UE_to_sector(u_,j) = another_interfering_pathloss_UE_to_sector(u_,j) + another_interfering_penetrationloss_UE_to_sector(u_,j) + another_interfering_shadow_fading_UE_to_sector(u_,j)-rx_antenna_gain(u_,j)-tx_antenna_gain(u_,j);
                    another_interfering_CL_UE_to_sector_linear(u_,j) = 10.^(0.1*another_interfering_CL_UE_to_sector(u_,j));
                    another_UE_up(u_,j) = element_gain_UE(u_,j)+rx_antenna_gain(u_,j);

                    another_BS_up(u_,j) =  element_gain_interference(u_,j);%+tx_antenna_gain(u_,j);



                    % ÿ��UE������书��
                    TX_power_ue = config.UE_tx_power;
                    %���ʿ���
                    P_min = 10^(-4)/1000;
                    R_min = P_min/TX_power_ue;
                    CL_x_ile = 88 + 10*log10(200/(config.bandwidth/1e6));
                    CL_x_ile_linear = 10.^(0.1*CL_x_ile);
                    Gama = config.Gama;
                    %�õ�����UE���غ�ķ���
                    another_interfering_TX_power(u_,j) = TX_power_ue * min(1,max(R_min,(another_interfering_CL_UE_to_sector_linear(u_,j)/CL_x_ile_linear).^Gama));
                end
            end

            obj.ul_ad_TX_power = another_interfering_TX_power;
            obj.ul_dl_UE_gain = another_UE_up;
            obj.ul_dl_BS_gain = another_BS_up;
            obj.ue_ue_pathloss = pl+penetration;
            if config.isACIR_loop %ѭ����¼ACIR
                loop_ACIR_dB = config.ACIR_lower:config.loop_step:config.ACIR_upper;
                loop_ACIR_linear = 10.^(0.1*loop_ACIR_dB);
                for i=1:length(loop_ACIR_dB)
                    loop_another_interfering_power(i,:) = sum(another_interfering_TX_power./another_interfering_CL_linear./ loop_ACIR_linear(i),2)';%�õ����յ���ϵͳ����UE��ĸ��Ź���֮��
                    ul_down_wideband_loop_SINR(i,:) = 10*log10(RX_power/(sum(loop_another_interfering_power(i,:))+sum(interfering_power(:))+thermal_noise_watts));%�������и������е�SINR
                end
                ul_down_wideband_loop_SINR(length(loop_ACIR_dB)+1,:) = 10*log10(RX_power/(sum(interfering_power(:))+thermal_noise_watts));%�����ACIR�����ʱ��SINR
            else
                ul_down_wideband_loop_SINR = 0;
            end
            ACIR_dB = config.ACIR_dB;%����У׼
            ACIR_linear = 10^(0.1*ACIR_dB);
            another_interfering_RX_power = sum(another_interfering_TX_power./another_interfering_CL_linear./ACIR_linear,2)'; %�õ����յ���ϵͳ����UE��ĸ��Ź���֮��
            % Interfering_power_asy = interfering_power(:)+another_interfering_RX_power(:);
            total_interfering_power=sum(interfering_power(:))+sum(another_interfering_RX_power(:));
            obj.ul_dl_wideband_SINR=10*log10(RX_power/(total_interfering_power+thermal_noise_watts));

        end

        %%  ���и������У�BS����BS��ֻ��˫ϵͳʱ�Ż�ʹ��
        function [user_macroscopic_pathloss,down_ul_wideband_loop_SINR] = down_up_model(obj,config,UEs,bs_beam_gain,ue_beam_gain,networkPathlossMap)
            interfering_eNodeBs = obj.attached_eNodeB.in_interf_eNodeB_sectors; %�õ�����eNodeB
            obj.penetration_loss = obj.downlink_channel.macroscopic_penetration_loss; %�õ���UE�������BS�Ĵ���
            user_macroscopic_pathloss        = obj.downlink_channel.macroscopic_pathloss;%_900;   %��+��Ԫ����
            user_shadow_fading_loss          = obj.downlink_channel.shadow_fading_pathloss;%�õ���UE�������BS����Ӱ˥��
            user_CL = user_macroscopic_pathloss + user_shadow_fading_loss + obj.penetration_loss -bs_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id)- ue_beam_gain(obj.id,obj.attached_eNodeB.eNodeB_id);%�����CL
            user_CL_linear = 10^(0.1*user_CL);%ת��Ϊ����ֵ

            %BS_antenna_gain = obj.downlink_channel.macroscopic_pathloss_model.BS_antenna_gain;
            % BS_antenna_gain = networkPathlossMap.BS_antenna_gain;
            % UE_antenna_gain = networkPathlossMap.UE_antenna_gain;

            the_RB_grid  = obj.downlink_channel.RB_grid;
            nRB          = floor(the_RB_grid.n_RB/config.n_UE_served_per_BS);
            nSC          = nRB*2;

            % ����
            thermal_noise_watts_per_half_RB = obj.thermal_noise_W_RB/2;
            thermal_noise_watts = thermal_noise_watts_per_half_RB*nSC;

            interfering_ue_id = [];
            tmp2 = [obj.attached_eNodeB.attached_UEs_vector.id]; %�õ���UE����BS�ķ����û���
            index_in_cell = find(tmp2 == obj.id);% index_in_cellΪUE��������С����UE������

            for s_ = 1:length(interfering_eNodeBs) %�õ����Ż�վ���ڷ�����û���
                tmp = [interfering_eNodeBs(s_).attached_UEs_vector.id];%�õ����Ż�վ�ķ����û���
                if length(tmp) == length(tmp2)
                    if isempty(interfering_ue_id) %�õ�����UE
                        interfering_ue_id = tmp(index_in_cell);
                    else
                        interfering_ue_id = [interfering_ue_id tmp(index_in_cell)];
                    end
                else

                    if isempty(interfering_ue_id)
                        interfering_ue_id = tmp(NR_randperm(length(tmp),1));
                    else
                        interfering_ue_id = [interfering_ue_id tmp(NR_randperm(length(tmp),1))];
                    end
                end
            end
            interfering_ue_id = sort(interfering_ue_id); % Ϊ����UE����

            UE_pos_vector = obj.downlink_channel.macroscopic_pathloss_model.UE_pos_vector;
            interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);% ����UE��λ��
            % ����UE����ǰUE�ķ���BS�������ģ���ϵͳ�����и��ţ�
            interfering_macroscopic_pathloss_UE = obj.downlink_channel.ul_interfering_macroscopic_pathloss(obj.attached_eNodeB.eNodeB_id,interfering_ue_id);%��+��Ԫ
            interfering_penetration_loss = obj.downlink_channel.ul_interfering_penetration_loss(obj.attached_eNodeB.eNodeB_id,interfering_ue_id);%��
            interfering_shadow_fading_loss = obj.downlink_channel.ul_interfering_shadow_fading_pathloss(obj.attached_eNodeB.parent_eNodeB.id,interfering_ue_id);%��˥
            interfering_CL = interfering_macroscopic_pathloss_UE + interfering_penetration_loss + interfering_shadow_fading_loss-bs_beam_gain(interfering_ue_id,obj.attached_eNodeB.eNodeB_id)'-ue_beam_gain(interfering_ue_id,obj.attached_eNodeB.eNodeB_id)';%������
            interfering_CL_linear = 10.^(0.1*interfering_CL);

            % ����UE�����Է���BS�������ģ����ڼ������UE�Ĺ��أ�
            interfering_sector_id = fix((interfering_ue_id-1)/config.UE_per_eNodeB)+1; %����������
            interfering_enb_id = obj.downlink_channel.macroscopic_pathloss_model.sector_idx_mapping(interfering_sector_id,1); %����eNodeB��
            interfering_enb_id = interfering_enb_id';

            for u_ = 1:length(interfering_ue_id)
                interfering_pathloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_pathloss_eNodeB(interfering_ue_pos(u_,:),interfering_sector_id(u_));%��+��Ԫ����

                interfering_UE_gain(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_ue_antenna_gain(interfering_ue_pos(u_,:),interfering_sector_id(u_));
                interfering_BS_gain(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_bs_antenna_gain(interfering_ue_pos(u_,:),interfering_sector_id(u_));
                interfering_penetrationloss_UE_to_sector(u_) = obj.downlink_channel.macroscopic_pathloss_model.get_penetration_loss_eNodeB(interfering_ue_pos(u_,:),interfering_sector_id(u_));%��͸���
                interfering_shadow_fading_UE_to_sector(u_) = obj.downlink_channel.shadow_fading_model.get_pathloss(interfering_ue_pos(u_,:),interfering_enb_id(u_));%��Ӱ˥��
                rx_antenna_gain(u_)=ue_beam_gain(interfering_ue_id(u_),UEs(interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);%UE�˲�����������
                beam_antenna_gain(u_)=bs_beam_gain(interfering_ue_id(u_),UEs(interfering_ue_id(u_)).attached_eNodeB.eNodeB_id);%BS�˲�����������

                interfering_gain_UE(u_) = interfering_UE_gain(u_); %+ rx_antenna_gain(u_);
                interfering_gain_BS(u_)  = interfering_BS_gain(u_) + beam_antenna_gain(u_);
            end
            interfering_CL_UE_to_sector = interfering_pathloss_UE_to_sector + interfering_penetrationloss_UE_to_sector + interfering_shadow_fading_UE_to_sector-rx_antenna_gain-beam_antenna_gain;%����UE�������BS��������
            interfering_CL_UE_to_sector_linear = 10.^(0.1*interfering_CL_UE_to_sector);

            obj.interfering_BS = interfering_gain_BS;
            obj.interfering_UE = interfering_gain_UE;

            % ÿ��UE����书��
            TX_power_ue = config.UE_tx_power;
            % ���ʿ���
            P_min = 10^(-4)/1000;
            R_min = P_min/TX_power_ue;
            CL_x_ile = 88 + 10*log10(200/(config.bandwidth/1e6));
            CL_x_ile_linear = 10.^(0.1*CL_x_ile);
            Gama = config.Gama;
            TX_power_PC = TX_power_ue * min(1,max(R_min,(user_CL_linear/CL_x_ile_linear)^Gama));%��UE�ķ��书��

            % Ϊ��UE���빦�����
            TX_power_PC = TX_power_PC * 10.^(0.1*obj.PCe);
            interfering_TX_power_PC = TX_power_ue .* min(1,max(R_min,(interfering_CL_UE_to_sector_linear./CL_x_ile_linear).^Gama)); %����UE�ķ��书��
            for u_=1:length(interfering_TX_power_PC)
                interfering_TX_power_PC(u_) = interfering_TX_power_PC(u_) * 10.^(0.1*UEs(interfering_ue_id(u_)).PCe);%Ϊ����UE���빦�����
            end

            obj.UE_tx_power = TX_power_PC;
            RX_power = TX_power_PC/user_CL_linear; %�����UE����BS�Ľ��չ���
            interfering_RX_power = interfering_TX_power_PC./interfering_CL_linear; %������յ�ͬһϵͳ����������UE�Ĺ���

            another_interfering_eNodeBs = obj.attached_eNodeB.ad_interf_eNodeB_sectors;%�õ���ϵͳ��������
            another_interfering_eNodeB_ids = [another_interfering_eNodeBs.eNodeB_id]; %�ҵ��ڶ���ϵͳ�ĸ�������id��
            for enb=1:length(another_interfering_eNodeB_ids)
                %����ÿ������BS�۲⵽����BS��ˮƽ�ǣ���ת��Ϊ(-180,180]
                phi_down_up(enb)=atan2(another_interfering_eNodeBs(enb).parent_eNodeB.pos(2)-obj.attached_eNodeB.parent_eNodeB.pos(2),another_interfering_eNodeBs(enb).parent_eNodeB.pos(1)-obj.attached_eNodeB.parent_eNodeB.pos(1))./pi*180-obj.attached_eNodeB.azimuth;
                phi_down_up(enb)=utils.miscUtils.wrapToAll180(phi_down_up(enb));
                %����ÿ������UE����UE�Ĵ�ֱ��
                switch config.scene_type
                    case {'RMa','UMA','UMI','InH','InH2','InH3'} %ͬ����˫ϵͳվַ�߶���ͬ
                        h1(enb)=config.site_height; %h1Ϊvictim��վ��
                        h2(enb)=config.site_height; % h2Ϊaggressorϵͳ��վ��
                        another_tx_power(enb)=config.eNodeB_tx_power/config.n_UE_served_per_BS;%���书��������йأ��̶���ͬʱ������û����й�
                    otherwise %�칹��վַ�߶Ȳ�ͬ
                        if another_interfering_eNodeB_ids(enb)> networkPathlossMap.num_first_sectors %���Ż�վ�ǵڶ���ϵͳ��(�߶Ȳ�Ϊ�ڶ���ϵͳվ���ڵ�һ��ϵͳ��վ�߲�ֵ)
                            h1(enb)=config.site_height;
                            h2(enb)=config.site_height2;
                            another_tx_power(enb)=config.eNodeB_tx_power2/config.n_UE_served_per_BS;
                        else %���Ż�վ�ǵ�һ��ϵͳ��(�߶Ȳ�Ϊ��һ��ϵͳվ���ڵڶ���ϵͳ��վ�߲�ֵ)
                            h1(enb)=config.site_height2;
                            h2(enb)=config.site_height;
                            another_tx_power(enb)=config.eNodeB_tx_power/config.n_UE_served_per_BS;
                        end
                end
                %����BS��۲⵽����BS�Ĵ�ֱ��
                theta_down_up(enb) = atan2(sqrt((obj.attached_eNodeB.parent_eNodeB.pos(1)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(1)).^2+(obj.attached_eNodeB.parent_eNodeB.pos(2)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(2))^2), h2(enb)- h1(enb))./pi*180;
                %����ÿ������BS��Ը�UE��������BS�ĵ�Ԫ����
                element_gain_interference(enb) = another_interfering_eNodeBs(enb).antenna.elementPattern(theta_down_up(enb),phi_down_up(enb));
                %�����UE��������BS�����ÿ������BS�ĵ�Ԫ����
                theta_opposite_observation(enb)=atan2(sqrt((obj.attached_eNodeB.parent_eNodeB.pos(1)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(1)).^2+(obj.attached_eNodeB.parent_eNodeB.pos(2)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(2))^2), h1(enb)- h2(enb))./pi*180;
                phi_opposite_observation(enb)=atan2(obj.attached_eNodeB.parent_eNodeB.pos(2)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(2),obj.attached_eNodeB.parent_eNodeB.pos(1)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(1))./pi*180-another_interfering_eNodeBs(enb).azimuth;

                %% �Ƕȵ���Ϊ(-180,180]
                phi_opposite_observation(enb)=utils.miscUtils.wrapToAll180(phi_opposite_observation(enb));
                %����victimϵͳ��վ����ڸ���BS�ĵ�Ԫ����
                element_gain_victim(enb) = obj.attached_eNodeB.antenna.elementPattern(theta_opposite_observation(enb),phi_opposite_observation(enb));

                %������ߵ�3d����
                d_3d(enb) = sqrt((obj.attached_eNodeB.parent_eNodeB.pos(1)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(1)).^2+(obj.attached_eNodeB.parent_eNodeB.pos(2)-another_interfering_eNodeBs(enb).parent_eNodeB.pos(2)).^2+(h2(enb)-h1(enb)).^2);
                if obj.id <= networkPathlossMap.num_first_UEs
                    frequency=config.frequency2;
                else
                    frequency=config.frequency;
                end
                if strcmp(config.macroscopic_pathloss_model,'TS38900') && strcmp(config.macroscopic_pathloss_model2,'TS38900') %�������ϵͳ����NR������NR�ļ��㹫ʽ
                    switch config.scene_type %������NLOS·��ģ�ͣ��칹ʱ���ô󳡾��Ĵ�����Ĺ�ʽ
                        case {'UMA','UMa_to_InH','UMa_to_UMi'}
                            temp_model = macroscopic_pathloss_models.TS38900PathlossModel(frequency,'urban_macro',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                UMA_std = 6;%UMA NLOS ����Ӱ˥���׼��
                            else
                                UMA_std = 4;%UMA LOS ����Ӱ˥���׼��
                            end
                            sf_UMA(enb) = normrnd(0,UMA_std); %������̬�ֲ�
                            pl(enb) = pl(enb) + sf_UMA(enb);
                            %                               pl(enb)=13.54+39.08*log10(d_3d(enb))+20*log10(frequency/1e9)-0.6*(obj.height-1.5)+sf_UMA(enb);
                        case {'UMI','UMi_to_InH'}
                            temp_model = macroscopic_pathloss_models.TS38900PathlossModel(frequency,'urban_micro',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                UMI_std = 7.82;%UMI NLOS ����Ӱ˥���׼��
                            else
                                UMI_std = 4;%UMI LOS ����Ӱ˥���׼��
                            end
                            sf_UMI(enb) = normrnd(0,UMI_std); %������̬�ֲ�
                            pl(enb) = pl(enb) + sf_UMI(enb);
                            %                               pl(enb)=35.3*log10(d_3d(enb))+22.4+21.3*log10(frequency/1e9)-0.3*(obj.height-1.5)+sf_UMI(enb);
                        case {'InH','InH2','InH3'}
                            temp_model = macroscopic_pathloss_models.TS38900PathlossModel(frequency,'indoor',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                InH_std = 8.03;%UMA NLOS ����Ӱ˥���׼��
                            else
                                InH_std = 3;%UMA NLOS ����Ӱ˥���׼��
                            end
                            sf_InH(enb) = normrnd(0,InH_std); %������̬�ֲ�
                            pl(enb) = pl(enb) + sf_InH(enb);
                            %                               pl(enb)=17.3+38.3*log10(d_3d(enb))+24.9*log10(frequency/1e9)+sf_InH(enb);%NLOS
                            %                           pl(enb)=32.4+17.3*log10(d_3d(enb))+20*log10(frequency/1e9);%LOS
                        case 'RMa'
                            temp_model = macroscopic_pathloss_models.TS38900PathlossModel(frequency,'rural_macro',config);
                            temp_model.h_ut = 35;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                RMA_std = 8;%UMA NLOS ����Ӱ˥���׼��
                            else
                                RMA_std = 4;%UMA NLOS ����Ӱ˥���׼��
                            end
                            sf_RMA(enb) = normrnd(0,RMA_std); %������̬�ֲ�
                            pl(enb) = pl(enb) + sf_RMA(enb);
                            %                               pl(enb)=161.04-7.1*log10(20)+7.5*log10(5)-(24.37-3.7*(5/35).^2)*log10(35)+(43.42-3.1*log10(35))*log10(d_3d(enb)-3)+20*log10(frequency/1e9)-(3.2*(log10(11.75*obj.height)).^2-4.97)+sf_RMa(enb);
                    end
                elseif strcmp(config.macroscopic_pathloss_model,'TS38901') && strcmp(config.macroscopic_pathloss_model2,'TS38901')
                    switch config.scene_type %������NLOS·��ģ�ͣ��칹ʱ���ô󳡾��Ĵ�����Ĺ�ʽ
                        case {'UMA','UMa_to_InH','UMa_to_UMi'}
                            temp_model = macroscopic_pathloss_models.TS38901PathlossModel(frequency,'urban_macro',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                UMA_std = 6;%UMA NLOS ����Ӱ˥���׼��
                            else
                                UMA_std = 4;%UMA NLOS ����Ӱ˥���׼��
                            end
                            sf_UMA(enb) = normrnd(0,UMA_std); %������̬�ֲ�
                            pl(enb) = pl(enb) + sf_UMA(enb);
                            %                               pl(enb)=13.54+39.08*log10(d_3d(enb))+20*log10(frequency/1e9)-0.6*(obj.height-1.5)+sf_UMA(enb);
                        case {'UMI','UMi_to_InH'}
                            temp_model = macroscopic_pathloss_models.TS38901PathlossModel(frequency,'urban_micro',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                UMI_std = 7.82;%UMA NLOS ����Ӱ˥���׼��
                            else
                                UMI_std = 4;%UMA NLOS ����Ӱ˥���׼��
                            end
                            sf_UMI(enb) = normrnd(0,UMI_std); %������̬�ֲ�
                            pl(enb) = pl(enb) + sf_UMI(enb);
                            %                               pl(enb)=35.3*log10(d_3d(enb))+22.4+21.3*log10(frequency/1e9)-0.3*(obj.height-1.5)+sf_UMI(enb);
                        case {'InH','InH2','InH3'}
                            temp_model = macroscopic_pathloss_models.TS38901PathlossModel(frequency,'indoor',config);
                            temp_model.h_ut = 25;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                InH_std = 8.03;%UMA NLOS ����Ӱ˥���׼��
                            else
                                InH_std = 3;%UMA NLOS ����Ӱ˥���׼��
                            end
                            sf_InH(enb) = normrnd(0,InH_std); %������̬�ֲ�
                            pl(enb) = pl(enb) + sf_InH(enb);
                            %                               pl(enb)=17.3+38.3*log10(d_3d(enb))+24.9*log10(frequency/1e9)+sf_InH(enb);%NLOS
                            %                           pl(enb)=32.4+17.3*log10(d_3d(enb))+20*log10(frequency/1e9);%LOS
                        case 'RMa'
                            temp_model = macroscopic_pathloss_models.TS38901PathlossModel(frequency,'rural_macro',config);
                            temp_model.h_ut = 35;
                            [pl(enb) comp] = temp_model.pathloss(d_3d(enb),25,0);
                            if comp
                                RMA_std = 8;%UMA NLOS ����Ӱ˥���׼��
                            else
                                RMA_std = 4;%UMA NLOS ����Ӱ˥���׼��
                            end
                            sf_RMA(enb) = normrnd(0,RMA_std); %������̬�ֲ�
                            pl(enb) = pl(enb) + sf_RMA(enb);
                            %                               pl(enb)=161.04-7.1*log10(20)+7.5*log10(5)-(24.37-3.7*(5/35).^2)*log10(35)+(43.42-3.1*log10(35))*log10(d_3d(enb)-3)+20*log10(frequency/1e9)-(3.2*(log10(11.75*obj.height)).^2-4.97)+sf_RMa(enb);
                    end
                else %UMA��RMA����ΪLTE
                    LTE_std=10;%LTE����Ӱ˥���׼��Ϊ10dB
                    LTE_sf(enb) = normrnd(0,LTE_std); %������̬�ֲ�
                    switch config.scene_type
                        case 'RMa'%�����RMa����
                            pl(enb) = 69.55+26.16*log10(config.frequency/1e6)-13.82*log10(config.site_height)+(44.9-6.55*log10(config.site_height))*log10(d_3d(enb)/1e3)-4.78*(log10(config.frequency/1e6))^2+18.33*log10(config.frequency/1e6)-40.94+LTE_sf(enb);
                        otherwise %�����UMa�������߰���UMA���칹����
                            pl(enb) = 40*(1-4e-3*abs(h2(enb)-h1(enb)))*log10(d_3d(enb)/1e3)-18*log10(abs(h2(enb)-h1(enb)))+21*log10(config.frequency/1e6)+80+LTE_sf(enb);
                    end

                end
                switch config.scene_type
                    case {'UMA','UMI','RMa','InH','UMa_to_UMi','InH2'}
                        penetration(enb)=0;%��վ���Ż�վ�޴�͸���
                    case {'UMi_to_InH','UMa_to_InH'} %��������ڵ��칹�������Ǵ�͸һ��ǽ
                        L_glass = 2+0.2*config.frequency/1e9;
                        L_concrete = 5+4*config.frequency/1e9;
                        N_low = 0+4.4*randn(1);
                        low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_low;
                        L_IRRglass = 23+0.3*config.frequency/1e9;
                        N_high = 0+6.5*randn(1);
                        high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_high;
                        if rand(1)>0.5 % 50%�ʹ���50%�ߴ���
                            penetration(enb)=low;
                        else
                            penetration(enb)=high;
                        end
                    case {'InH3'} %��������ڵ��칹�������Ǵ�͸һ��ǽ
                        L_glass = 2+0.2*config.frequency/1e9;
                        L_concrete = 5+4*config.frequency/1e9;
                        N_low = 0+4.4*randn(1);
                        low = 5-10*log10(0.3*10^(-L_glass/10)+0.7*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_low;
                        L_IRRglass = 23+0.3*config.frequency/1e9;
                        N_high = 0+6.5*randn(1);
                        high = 5-10*log10(0.7*10^(-L_IRRglass/10)+0.3*10^(-L_concrete/10))+0.5*min(25*rand(1,1),25*rand(1,1))+N_high;
                        if rand(1)>0.5 % 50%�ʹ���50%�ߴ���
                            penetration(enb)=low*2;
                        else
                            penetration(enb)=high*2;
                        end
                end
                %������ϵͳ����BS����UE����BS��������
                another_interfering_CL(enb) = pl(enb) + penetration(enb)-element_gain_interference(enb) - element_gain_victim(enb) - ue_beam_gain(obj.id,another_interfering_eNodeB_ids(enb)) - bs_beam_gain(obj.id,another_interfering_eNodeB_ids(enb));
                another_interfering_CL_sent(enb) = pl(enb) + penetration(enb)-element_gain_victim(enb)  - ue_beam_gain(obj.id,another_interfering_eNodeB_ids(enb)) ;
                another_UE(enb) = -element_gain_interference(enb)-ue_beam_gain(obj.id,another_interfering_eNodeB_ids(enb)) ;
                another_BS(enb) =  element_gain_victim(enb)+bs_beam_gain(obj.id,another_interfering_eNodeB_ids(enb));
                another_interfering_CL_inter(enb) = pl(enb) + penetration(enb)-element_gain_interference(enb)  - bs_beam_gain(obj.id,another_interfering_eNodeB_ids(enb));
                another_interfering_CL_linear(enb) = 10.^(0.1*another_interfering_CL(enb)); %ת��Ϊ����ֵ
            end
            %               obj.bs_bs_pathloss = pl + penetration;%·��
            obj.dl_ul_UE_gain = another_UE;
            obj.dl_ul_BS_gain = another_BS;
            %               obj.bs_bs_pathloss = another_interfering_CL_sent;%·��+������������
            obj.bs_bs_pathloss = another_interfering_CL_inter ;%·��+������������

            if config.isACIR_loop %�Ƿ�ѭ��ͳ�Ʋ�ͬACIR�µĽ��
                loop_ACIR_dB = config.ACIR_lower:config.loop_step:config.ACIR_upper;
                loop_ACIR_linear = 10.^(0.1*loop_ACIR_dB);
                for i=1:length(loop_ACIR_dB)
                    loop_another_interfering_power(i,:) = another_tx_power./another_interfering_CL_linear./ loop_ACIR_linear(i);%�����յ�����ϵͳ����BSй¶�Ĺ���
                    %���Ż�վ�Թ�ַ����ϵͳ��վ�ĸ��Ű�·��ʽ�����޴�������蹲ַ������ϵͳ��վ��ȫ����
                    for j=1:length(loop_another_interfering_power(i,:))
                        if loop_another_interfering_power(i,j)==Inf || isnan(loop_another_interfering_power(i,j))
                            loop_another_interfering_power(i,j)=0; %��ϵͳ��ַBS�Ա�ϵͳ��BS�ĸ�������Ϊ0
                        end
                    end
                    down_ul_wideband_loop_SINR(i,:) = 10*log10(RX_power/(sum(loop_another_interfering_power(i,:))+sum(interfering_RX_power(:))+thermal_noise_watts));%�������и������е�SINR
                end
                down_ul_wideband_loop_SINR(length(loop_ACIR_dB)+1,:) = 10*log10(RX_power/(sum(interfering_RX_power(:))+thermal_noise_watts));
            else
                down_ul_wideband_loop_SINR = 0;
            end
            ACIR_dB = config.ACIR_dB; %У׼ʹ�õ�ACIRֵ
            ACIR_linear = 10^(0.1*ACIR_dB);
            another_interfering_RX_power=another_tx_power./another_interfering_CL_linear./ACIR_linear; %�õ����յ��ĸ���UE����;
            %���Ż�վ�Թ�ַ����ϵͳ��վ�ĸ��Ű�·��ʽ�����޴�������蹲ַ������ϵͳ��վ��ȫ����
            total_interfering_power=sum(interfering_RX_power(:))+sum(another_interfering_RX_power(find(another_interfering_RX_power~=inf & ~isnan(another_interfering_RX_power))));
            obj.dl_ul_wideband_SINR=10*log10(RX_power/(total_interfering_power+thermal_noise_watts));  %��У׼ʱʹ�õ�ACIR�µ�SINR��¼����
        end

        function dummy_up_down_model(obj)
            obj.ul_dl_wideband_SINR = NaN;
        end

        function dummy_down_up_model(obj)
            obj.dl_ul_wideband_SINR = NaN;
        end

        function struct_out = basic_information_in_struct(obj)
            struct_out.id                    = obj.id;
            struct_out.pos                   = obj.pos;
        end

        function distance_to_site = distance_to_attached_site(obj)
            % Return the distance to the attached site (function added for convenience)
            distance_to_site = sqrt(sum((obj.attached_site.pos -obj.pos).^2));
        end

        % Clear all non-basic info and leaves just basic information describing the UE
        function clear_non_basic_info(obj)
            obj.attached_site             = [];
            obj.attached_sector_idx       = [];
            obj.attached_eNodeB           = [];
            obj.walking_model             = [];
            obj.downlink_channel          = [];
            obj.RB_grid                   = [];
            obj.uplink_channel            = [];
            obj.clock                     = [];
            obj.wideband_SINR             = [];
            obj.ul_wideband_SINR          = [];
            obj.ul_dl_wideband_SINR       = [];
            obj.dl_ul_wideband_SINR       = [];
            obj.SINR_TTI                  = [];
            obj.ul_SINR_TTI               = [];
            obj.ul_dl_SINR_TTI            = [];
            obj.dl_ul_SINR_TTI            = [];
            obj.attached_eNodeB_id_TTI    = [];
            obj.deactivate_UE             = [];
            obj.receiver_noise_figure     = [];
            obj.thermal_noise_W_RB        = [];
            obj.penetration_loss          = [];
            obj.rx_antenna_gain           = [];
            obj.beam_antenna_gain         = [];
        end

    end
end
