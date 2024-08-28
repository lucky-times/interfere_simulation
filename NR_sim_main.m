function [output_results_file] = NR_sim_main(SYS_config)
%ϵͳ�������������
%
%���������
%SYS_config���������
%
%���������
%output_results_file���������ļ�

ticID_start_sim_begin = tic; % ����ʼ��ʱ��

%% ��ֹ�������ò���ȷ
SYS_config = NR_load_params_again( SYS_config );

%% ������Եȼ�
DEBUG_LEVEL = SYS_config.debug_level;

%% ������վ�����õ�·��ͼ������Ӱ˥��
[eNodeB_sites,eNodeB_sectors,networkPathlossMap,networkShadowFadingMap] = NR_generate_network(SYS_config);

%% ʱ��
networkClock = network_elements.clock(SYS_config.TTI_length);

for b_=1:length(eNodeB_sites)
    eNodeB_sites(b_).clock = networkClock;
end

%% �����û�
[UEs] = NR_generate_users(SYS_config,eNodeB_sites,eNodeB_sectors,networkPathlossMap,networkClock);
eNodeB_pos = [];
for i = 1:length(eNodeB_sites)
    eNodeB_pos = [eNodeB_pos;eNodeB_sites(i).pos];
end
UE_pos = [];
for i = 1:length(UEs)
    UE_pos = [UE_pos;UEs(i).pos];
end
figure;
scatter(eNodeB_pos(:,1), eNodeB_pos(:, 2));
hold on;
scatter(UE_pos(:,1), UE_pos(:, 2));
hold off;
axis equal;
title('���˷ֲ�ͼ');

if isempty(UEs)
    warning('No UEs generated.');
else  
     %% UE��������ĳ�ʼ�������ڼ�������
    for u_=1:length(UEs)
        
        % ���������ŵ� (����������ģ���͸��ĺ���Ӱ˥��)
        UEs(u_).downlink_channel = channel_models.downlinkChannelModel(UEs(u_));
                       
        % ���봫�����
        UEs(u_).downlink_channel.set_macroscopic_pathloss_model(networkPathlossMap);
        
        % ������Ӱ˥��
        if SYS_config.macroscopic_pathloss_is_model
            UEs(u_).downlink_channel.set_shadow_fading_model(networkShadowFadingMap);
        end
        
        % ��ȡUE����������RB������߶ȡ���ת�ǣ�
        UEs(u_).RB_grid = UEs(u_).downlink_channel.RB_grid;
%         UEs(u_).height = UEs(u_).downlink_channel.UE_height;
        UEs(u_).height = SYS_config.UE_height;
        UEs(u_).orientation = 0;
%         UEs(u_).orientation = UEs(u_).downlink_channel.UE_orientation;
    end
    
    %% ��ӡeNodeB
    if DEBUG_LEVEL>=2
        fprintf('eNodeB List\n');
        for b_=1:length(eNodeB_sites)
            eNodeB_sites(b_).print;
        end
        fprintf('\n');
    end
    
    %% ��ӡ�����û�
    if DEBUG_LEVEL>=2
        fprintf('User List\n');
        for u_=1:length(UEs)
            UEs(u_).print;
        end
        fprintf('\n');
    end
    
    %% ��ѭ��
    if DEBUG_LEVEL>=1
        fprintf('Entering main simulation loop, %5.0f TTIs\n',SYS_config.simulation_time_tti);
    end
    
    % ��ʼ��������
    ticID_start_sim_loop = tic;
    starting_time = toc(ticID_start_sim_loop);
    
    num_markers = 5;% �������ĳ���
    s_markings  = round(linspace(1,length(eNodeB_sectors),num_markers));%�������
    u_markings  = round(linspace(1,length(UEs),num_markers));%UE���
    
    % ʱ�ӳ�ʼ��
    while networkClock.current_TTI < SYS_config.simulation_time_tti
        networkClock.advance_1_TTI;
        if DEBUG_LEVEL>=1
            fprintf(' TTI %5.0f/%d: ',networkClock.current_TTI,SYS_config.simulation_time_tti);
        end
        
        % ��������
        if DEBUG_LEVEL>=1
            fprintf('\n               ');
        end
        if SYS_config.asynchronization %����ǲ�ͬ��
            if SYS_config.isDouble == true %�����˫ϵͳ�ҿ������첽����
                [ue_or_bs_gain_asy,ue_gain_asy,bs_or_ue_gain_asy,bs_gain_asy] = NR_beam_forming_asy(UEs,eNodeB_sectors,networkPathlossMap,SYS_config); %��ͬ��ʱ���ò�ͬ�Ĳ�������

                [beam_antenna_gain,beam_ue_gain] = NR_beam_forming(UEs,eNodeB_sectors,networkPathlossMap,SYS_config);

            else
                error('isDouble must be true')
            end
        else %�����ͬ������������Ĳ������ͺ���
            [beam_antenna_gain,beam_ue_gain] = NR_beam_forming(UEs,eNodeB_sectors,networkPathlossMap,SYS_config);
        end
       
        % ����SINR
        if DEBUG_LEVEL>=1
            fprintf('               ');
        end
        switch SYS_config.asynchronization_switch
          case 0
                %����
                kk_ = 1;
                for u_ = 1:length(UEs)
                    if ~UEs(u_).deactivate_UE
                        [signal_CL(kk_,:),interfering_CL(kk_,:),another_interfering_CL(kk_,:),UE_Rx_power(kk_),UE_Rx_intf_power(kk_),wideband_loop_SINR(:,:,kk_)] = UEs(u_).down_link_quality_model(SYS_config,beam_antenna_gain,beam_ue_gain, SYS_config);
                        kk_ = kk_ + 1;
                    else
                        UEs(u_).dummy_down_link_quality_model;
                    end
                    if ~isempty(find(u_==u_markings,1))
                        if DEBUG_LEVEL>=1
                            fprintf('+'); %������
                        end
                    end
                end
                for i = 1:length(UEs)
                    SINR(i) = UEs(i).wideband_SINR;
                    SINR2(i) = UEs(i).wideband_SINR2;
                    attached_sector_id(i) = UEs(i).attached_eNodeB.eNodeB_id;
                end
                figure;
                plot(attached_sector_id)
                x = linspace(0,360,length(UEs));
                ylabel('PCI')
                title('UAV����С����ţ����棩');
%                 [~,locs]= findpeaks(SINR);
%                 x_tick = x(locs);
%                 y_tick = SINR(locs);
                fig1 = figure;
                plot(x, SINR,LineWidth=2.5);
%                 hold on;
%                 plot(x_tick, y_tick, 'ro', LineWidth=2.5)
                xlabel('����/��')
                ylabel('SINR/dB')
%                 str = sprintf('������ά������㣬��վ��������%d��\nUAV������%dm��h=%dm', SYS_config.tilt, SYS_config.UE_r, SYS_config.UE_height);
                str = sprintf('������ά������㣬��վ���ߴ�ֱ����24��\nUAV������%dm��h=%dm', SYS_config.UE_r, SYS_config.UE_height);
                title(str)
% 
                str1 = sprintf('pos_r%d_h%d', SYS_config.UE_r, SYS_config.UE_height);
%                 saveas(fig1, str1, 'svg');        % ָ��·������

%                 fig2 = figure;
%                 plot(x, SINR2,LineWidth=2.5);
% %                 hold on;
% %                 plot(x_tick, y_tick, 'bo', LineWidth=2.5)
%                 xlabel('����/��')
%                 ylabel('SINR/dB')
% %                 str = sprintf('�������Ĺ�ʽ�������㣬��վ��������%d��\nUAV����=%dm��h=%dm', SYS_config.tilt, SYS_config.UE_r, SYS_config.UE_height);
%                 str = sprintf('�������Ĺ�ʽ�������㣬��վ���ߴ�ֱ����24��\nUAV����=%dm��h=%dm', SYS_config.UE_r, SYS_config.UE_height);
%                 title(str)
% 
%                 str2 = sprintf('alpha_r%d_h%d', SYS_config.UE_r, SYS_config.UE_height);
%                 saveas(fig2, str2, 'svg');        % ָ��·������


                error('��ͼ�ѽ���')

                
            case 1
                if SYS_config.isDouble == true
                    kk_ = 1;% ������¼activate_UE�������ݵ��±�
                    for u_ = 1:length(UEs)
                        if ~UEs(u_).deactivate_UE
                            [user_macroscopic_pathloss(kk_),ul_dl_wideband_loop_SINR(:,:,kk_)] = UEs(u_).up_down_model(SYS_config,UEs,ue_or_bs_gain_asy,ue_gain_asy,networkPathlossMap);
                            [dl_user_macroscopic_pathloss(kk_),dl_ul_wideband_loop_SINR(:,:,kk_)] = UEs(u_).down_up_model(SYS_config,UEs,bs_or_ue_gain_asy,bs_gain_asy,networkPathlossMap);
                            kk_ = kk_ + 1;
                        else
                            UEs(u_).dummy_up_down_model;
                            UEs(u_).dummy_down_up_model;
                        end
                        if ~isempty(find(u_==u_markings,1))
                            if DEBUG_LEVEL>=1
                                fprintf('*'); %������
                            end
                        end
                        % ����ÿ��TTI�õ���SINR
                        UEs(u_).save_SINR_in_TTI(SYS_config);
                    end
                end
               
            case 2
                UE_SINR_sum_dl=[];
                UE_SINR_sum_ul=[];
                %����
                kk_ = 1;
                for u_ = 1:length(UEs)
                    if ~UEs(u_).deactivate_UE
                        [signal_CL(kk_,:),interfering_CL(kk_,:),another_interfering_CL(kk_,:),UE_Rx_power(kk_),UE_Rx_intf_power(kk_),wideband_loop_SINR(:,:,kk_)] = UEs(u_).down_link_quality_model(SYS_config,beam_antenna_gain,beam_ue_gain);
                        kk_ = kk_ + 1;
                    else
                        UEs(u_).dummy_down_link_quality_model;
                    end
                    if ~isempty(find(u_==u_markings,1))
                        if DEBUG_LEVEL>=1
                            fprintf('+'); %������
                        end
                    end
                end

                % ����
                kk_ = 1;
                for u_ = 1:length(UEs)
                    if ~UEs(u_).deactivate_UE
                        [IoT(kk_),interference_level(kk_),BS_Rx_power(kk_),BS_Rx_intf_power(kk_),ul_wideband_loop_SINR(:,:,kk_)] = UEs(u_).up_link_quality_model(SYS_config,UEs,beam_antenna_gain,beam_ue_gain);
                        kk_ = kk_ + 1;
                    else
                        UEs(u_).dummy_up_link_quality_model;
                    end
                    if ~isempty(find(u_==u_markings,1))
                        if DEBUG_LEVEL>=1
                            fprintf('-');  %������
                        end
                    end
                    % ����ÿ��TTI�õ���SINR
                    %UEs(u_).save_SINR_in_TTI(SYS_config);
                end
                
                %�첽���ż���
                kk_ = 1;% ������¼activate_UE�������ݵ��±�
                    for u_ = 1:length(UEs)
                        if ~UEs(u_).deactivate_UE
                            ul_dl_wideband_loop_SINR(:,:,kk_) = UEs(u_).up_down_model(SYS_config,UEs,ue_or_bs_gain_asy,ue_gain_asy,networkPathlossMap);
                            dl_ul_wideband_loop_SINR(:,:,kk_) = UEs(u_).down_up_model(SYS_config,UEs,bs_or_ue_gain_asy,bs_gain_asy,networkPathlossMap);
                            kk_ = kk_ + 1;
                        else
                            UEs(u_).dummy_up_down_model;
                            UEs(u_).dummy_down_up_model;
                        end
                    end
        end
        
%         if SYS_config.asynchronization %�����в�ͬ��
%             if SYS_config.isDouble == true
%                 kk_ = 1;% ������¼activate_UE�������ݵ��±�
%                 for u_ = 1:length(UEs)
%                     if ~UEs(u_).deactivate_UE
%                         ul_dl_wideband_loop_SINR(:,:,kk_) = UEs(u_).up_down_model(SYS_config,UEs,ue_or_bs_gain_asy,ue_gain_asy,networkPathlossMap);
%                         dl_ul_wideband_loop_SINR(:,:,kk_) = UEs(u_).down_up_model(SYS_config,UEs,bs_or_ue_gain_asy,bs_gain_asy,networkPathlossMap);
%                         kk_ = kk_ + 1;
%                     else
%                         UEs(u_).dummy_up_down_model;
%                         UEs(u_).dummy_down_up_model;
%                     end
%                     if ~isempty(find(u_==u_markings,1))
%                         if DEBUG_LEVEL>=1
%                             fprintf('*'); %������
%                         end
%                     end
%                     % ����ÿ��TTI�õ���SINR
%                     UEs(u_).save_SINR_in_TTI(SYS_config);
%                 end
%             end
%         else %������ͬ����������и������к����и������е�SINR
%             % ����
%             kk_ = 1;
%             for u_ = 1:length(UEs)
%                 if ~UEs(u_).deactivate_UE
%                     [signal_CL(kk_,:),interfering_CL(kk_,:),another_interfering_CL(kk_,:),UE_Rx_power(kk_),UE_Rx_intf_power(kk_),wideband_loop_SINR(:,:,kk_)] = UEs(u_).down_link_quality_model(SYS_config,beam_antenna_gain,beam_ue_gain);
%                     kk_ = kk_ + 1;
%                 else
%                     UEs(u_).dummy_down_link_quality_model;
%                 end
%                 if ~isempty(find(u_==u_markings,1))
%                     if DEBUG_LEVEL>=1
%                         fprintf('+'); %������
%                     end
%                 end
%             end
%             
%             % ����
%             kk_ = 1;
%             for u_ = 1:length(UEs)
%                 if ~UEs(u_).deactivate_UE
%                     [IoT(kk_),interference_level(kk_),BS_Rx_power(kk_),BS_Rx_intf_power(kk_),ul_wideband_loop_SINR(:,:,kk_)] = UEs(u_).up_link_quality_model(SYS_config,UEs,beam_antenna_gain,beam_ue_gain);
%                     kk_ = kk_ + 1;
%                 else
%                     UEs(u_).dummy_up_link_quality_model;
%                 end
%                 if ~isempty(find(u_==u_markings,1))
%                     if DEBUG_LEVEL>=1
%                         fprintf('-');  %������
%                     end
%                 end
%                 % ����ÿ��TTI�õ���SINR
%                 UEs(u_).save_SINR_in_TTI(SYS_config);
%             end
%         end
        
        if DEBUG_LEVEL>=1
            fprintf('\n');
        end
        
        if networkClock.current_TTI < SYS_config.simulation_time_tti
            % Move UEs
            if ~SYS_config.keep_UEs_still
                move_all_UEs(SYS_config,UEs,networkPathlossMap,eNodeB_sectors);
            end
        end
        
        if mod(networkClock.current_TTI,10)==0
            elapsed_time = toc(ticID_start_sim_loop);
            time_per_iteration = elapsed_time / networkClock.current_TTI;
            estimated_time_to_finish = (SYS_config.simulation_time_tti - networkClock.current_TTI)*time_per_iteration;
            estimated_time_to_finish_h = floor(estimated_time_to_finish/3600);
            estimated_time_to_finish_m = estimated_time_to_finish/60 - estimated_time_to_finish_h*60;
            fprintf('Time to finish: %3.0f hours and %3.2f minutes\n',estimated_time_to_finish_h,estimated_time_to_finish_m);
        end
    end
    
    %% ͳ��
    % ��ͳ�Ƶ�UE����deactivate_UE = false��
    if DEBUG_LEVEL>=1
        fprintf('Generating statistics\n');
    end
    if SYS_config.keep_UEs_still
        activate_UE_id = NR_get_activate_UE( UEs );
        switch SYS_config.asynchronization_switch
            case 0
               NR_plot_SINR_CDF(SYS_config,UEs,activate_UE_id); 
               NR_throughput(SYS_config,wideband_loop_SINR,ul_wideband_loop_SINR,activate_UE_id);
            case 1
               NR_plot_SINR_CDF_asy(SYS_config,UEs,activate_UE_id);
               NR_throughput_asy(SYS_config,ul_dl_wideband_loop_SINR,dl_ul_wideband_loop_SINR,activate_UE_id);
               NR_user_PL_gain(SYS_config,user_macroscopic_pathloss,dl_user_macroscopic_pathloss);
               NR_interfering_gain(SYS_config,UEs,activate_UE_id);
               NR_plot_UETxPower_CDF(SYS_config,UEs,activate_UE_id);
               %NR_user_dl_to_ul_PL_gain(SYS_config,dl_user_macroscopic_pathloss);
            case 2
               %ͬ������SINR
               all_SINR = [UEs.wideband_SINR];
               UE_SINR = all_SINR(activate_UE_id);
               UE_SINR_sum_dl=[UE_SINR_sum_dl UE_SINR];
               ul_all_SINR = [UEs.ul_wideband_SINR];
               UE_SINR_ul = ul_all_SINR(activate_UE_id);
               UE_SINR_sum_ul=[UE_SINR_sum_ul UE_SINR_ul];
               %�첽����SINR
               ul_dl_all_SINR = [UEs.ul_dl_wideband_SINR];
               UE_SINR = ul_dl_all_SINR(activate_UE_id);
               UE_SINR_sum_dl=[UE_SINR_sum_dl UE_SINR];
               dl_ul_all_SINR = [UEs.dl_ul_wideband_SINR];
               UE_SINR_ul = dl_ul_all_SINR(activate_UE_id);
               UE_SINR_sum_ul=[UE_SINR_sum_ul UE_SINR_ul];
               NR_plot_SINR_CDF_half(UE_SINR_sum_dl,UE_SINR_sum_ul);
               
               wideband_loop_SINR_sum=wideband_loop_SINR;
               ul_wideband_loop_SINR_sum=ul_wideband_loop_SINR;
               wideband_loop_SINR_sum=cat(3,wideband_loop_SINR_sum,ul_dl_wideband_loop_SINR);
               ul_wideband_loop_SINR_sum=cat(3,ul_wideband_loop_SINR_sum,dl_ul_wideband_loop_SINR);
               NR_throughput_half(SYS_config,wideband_loop_SINR_sum,ul_wideband_loop_SINR_sum,activate_UE_id);
        end
        if SYS_config.asynchronization
            NR_plot_bs_bs_pathloss(SYS_config,UEs,activate_UE_id);
            NR_plot_ue_ue_pathloss(SYS_config,UEs,activate_UE_id);
        end
        NR_plot_couplingloss_CDF(SYS_config,UEs,beam_antenna_gain,beam_ue_gain,activate_UE_id);
        


%
       % NR_plot_receive_power_CDF( SYS_config,ul_to_dl_Agg_BS_Tx,  ul_to_dl_Agg_Rx,dl_to_ul_Agg_BS_Tx,dl_to_ul_Agg_Rx);
        %NR_dl_ul_gain(SYS_config,UEs,activate_UE_id);
        %NR_dl_ul_gain(SYS_config,ue_or_bs_gain_asy);
        %NR_throughput_asy(SYS_config,ul_dl_wideband_loop_SINR,dl_ul_wideband_loop_SINR,activate_UE_id);
      %  NR_plot_signal_and_intf_CL_CDF(,interfering_CL,another_interfering_CL,activate_UE_id);
       % NR_plot_receive_power_CDF( SYS_config,UE_Rx_power,UE_Rx_intf_power,BS_Rx_power,BS_Rx_intf_power );
        %NR_plot_receive_power_CDF(SYS_config,BS_Tx_power_asy,UE_RX_power_asy,BS_Tx_intf_power_asy,UE_Rx_intf_power_asy)
    %    NR_plot_UETxPower_CDF(SYS_config,UEs,activate_UE_id);
%             %����IoTcdf
%         if ~SYS_config.asynchronization
%             %������������
%             NR_plot_couplingloss_CDF(SYS_config,UEs,beam_antenna_gain,beam_ue_gain,activate_UE_id);
%             %���湦��cdf���
%             NR_plot_UETxPower_CDF(SYS_config,UEs,activate_UE_id);
%             %������ϵͳ������cdf���
%             NR_plot_signal_and_intf_CL_CDF(signal_CL,interfering_CL,another_interfering_CL,activate_UE_id);
%             %����IoTcdf
%             NR_plot_IoT_CDF(SYS_config,IoT);
%             % ����blocking cdf
%             NR_plot_blocking_CDF( SYS_config,interference_level );
%             % ���ƽ��յ㹦�ʼ���������
%             NR_plot_receive_power_CDF( SYS_config,UE_Rx_power,UE_Rx_intf_power,BS_Rx_power,BS_Rx_intf_power );
%             % ��ͬACIR�µ��������仯
%             if SYS_config.isACIR_loop && SYS_config.isDouble
%                 NR_throughput(SYS_config,wideband_loop_SINR,ul_wideband_loop_SINR,activate_UE_id);
%             end
%             % �̶�ACIR�µ�����������ӡ��Command Windows��
%             NR_plot_thput_CDF(SYS_config,UEs,activate_UE_id);
%             % ����SINR��CDF
%             NR_plot_SINR_CDF(SYS_config,UEs,activate_UE_id);
%         else
%             % �첽ʱ��
%             NR_plot_thput_CDF_asy(SYS_config,UEs,activate_UE_id);
%             if SYS_config.isACIR_loop && SYS_config.isDouble
%                 NR_throughput_asy(SYS_config,ul_dl_wideband_loop_SINR,dl_ul_wideband_loop_SINR,activate_UE_id);
%             end
%         end
    end
    
    %% ������ʱ
    finish_time_s_full = toc(ticID_start_sim_begin);
    finish_time_m = floor(finish_time_s_full/60);
    finish_time_s = finish_time_s_full-finish_time_m*60;
    if DEBUG_LEVEL>=1
        fprintf('Simulation finished\n');
        fprintf(' Total elapsed time: %.0fm, %.0fs\n',finish_time_m,finish_time_s);
    end
    
    if DEBUG_LEVEL>=1
        fprintf('Saving results to %s\n',SYS_config.results_file);
    end
    
    %% ����ʱ�����ĸ������ݵı���
    % Compact (or not) results
    if SYS_config.compact_results_file
        
        networkPathlossMap.delete_everything_except_cell_assignments();
        
        clear the_UE_traces;

        % ��UE�Ĳ������Ա�ɽṹ�屣��
        for u_=1:length(UEs)
            the_UE_traces(u_).position        = UEs(u_).pos;
            the_UE_traces(u_).attached_site   = UEs(u_).attached_site.id;
            the_UE_traces(u_).attached_eNodeB = UEs(u_).attached_eNodeB.eNodeB_id;
            the_UE_traces(u_).attached_eNodeB_id_TTI = UEs(u_).attached_eNodeB_id_TTI;
            the_UE_traces(u_).SINR_TTI = UEs(u_).SINR_TTI;
            the_UE_traces(u_).ul_SINR_TTI = UEs(u_).ul_SINR_TTI;
            the_UE_traces(u_).ul_dl_SINR_TTI = UEs(u_).ul_dl_SINR_TTI;
            the_UE_traces(u_).dl_ul_SINR_TTI = UEs(u_).dl_ul_SINR_TTI;
            
            UEs_struct(u_)                    = UEs(u_).basic_information_in_struct();
            UEs(u_).clear_non_basic_info();
        end
        UEs = UEs_struct; 
        % ������������
        %save(fullfile(SYS_config.results_folder,SYS_config.results_file),'SYS_config','networkPathlossMap','eNodeB_sites','eNodeB_sectors','UEs','the_UE_traces','finish_time_s_full','-v7.3');
    else
        
        if SYS_config.delete_pathloss_at_end
            networkPathlossMap.pathloss = [];
            if exist('networkShadowFadingMap','var')
                networkPathlossMap.pathloss = [];
            else
                % Do nothing
            end
        end
        save(fullfile(SYS_config.results_folder,SYS_config.results_file),'-v7.3');
        if DEBUG_LEVEL>=1
            fprintf('Traces saved in the standard format\n');
        end
    end
end

output_results_file = fullfile(SYS_config.results_folder,SYS_config.results_file);

utils.miscUtils.tidy_up_memory_before_closing(UEs,eNodeB_sectors,eNodeB_sites);
clear the_UE_traces
clear networkPathlossMap
clear networkShadowFadingMap
clear UEs
clear eNodeBs
clear eNodeBs_sectors

function move_all_UEs(SYS_config,UEs,networkPathlossMap,eNodeBs_sectors)
%�ú���������չUE�ƶ�����

some_UE_out_of_ROI_this_TTI = false;
num_first_UEs = networkPathlossMap.num_first_UEs;
num_first_sectors = networkPathlossMap.num_first_sectors;
for u_ = 1:length(UEs)
    UEs(u_).move;
    [ x_range,y_range ] = networkPathlossMap.valid_range;
    
    % ���UE�߳���ROI������������ƶ���ROI������λ��
    %
    % ��ͬʱҲʵ����һ���򵥵��л����򣬿����ڴ˻����Ͻ�������
    
    ROI_teleport       = ~UEs(u_).is_in_roi(SYS_config,x_range,y_range);% UE�Ƿ���ROI�У���Ļ�ROI_teleport=true
    handover_requested = UEs(u_).cell_change.requested;

    if ROI_teleport || handover_requested
        old_eNodeB_id = UEs(u_).attached_eNodeB.eNodeB_id;
        if ROI_teleport %���������ROI
            new_UE_position = networkPathlossMap.random_position;%���ѡ��λ��
            
            [new_site_id,new_sector_id,new_eNodeB_id] = networkPathlossMap.cell_assignment(SYS_config,u_,new_UE_position,num_first_UEs);
            
            %λ�ø���
            UEs(u_).pos = new_UE_position;
        elseif handover_requested
            new_eNodeB_id                     = UEs(u_).cell_change.target_eNodeB;
            UEs(u_).cell_change.requested     = false; % �����л���
            UEs(u_).cell_change.target_eNodeB = [];
        end
        
         % ��UE��ԭ��վ���룬���������»�վ
        UEs(u_).start_handover(eNodeBs_sectors(new_eNodeB_id));
        

        % ͳ�ƻ�ͳ��UE��SINR
        if SYS_config.isWraparound
            compute_only_UEs_from_this_eNodeBs = 1:57;
        else
            if SYS_config.isDouble
                compute_only_UEs_from_this_eNodeBs = 1:num_first_sectors;
            else
                compute_only_UEs_from_this_eNodeBs = 1:length(eNodeBs_sectors);
            end
        end
        if isempty(find(UEs(u_).attached_eNodeB.eNodeB_id == compute_only_UEs_from_this_eNodeBs,1))
            % ��ͳ�Ƹ�UE
            UEs(u_).deactivate_UE = true;
        else
            % ͳ�Ƹ�UE
            UEs(u_).deactivate_UE = false;
        end
        
        % ��ӡһЩ������Ϣ
        if ~some_UE_out_of_ROI_this_TTI
            if SYS_config.debug_level>=1
                fprintf(1,'\n');
            end
            some_UE_out_of_ROI_this_TTI = true;
        end
        if SYS_config.debug_level>=1
            if ROI_teleport
                fprintf('              UE %g going out of ROI, teleporting to %g %g. eNodeB %g -> eNodeB %g\n',UEs(u_).id,new_UE_position(1),new_UE_position(2),old_eNodeB_id,new_eNodeB_id);
            elseif handover_requested
                fprintf('              UE %g handover request. eNodeB %g -> eNodeB %g\n',UEs(u_).id,old_eNodeB_id,new_eNodeB_id);
            end
        end
    end
end
if some_UE_out_of_ROI_this_TTI
    if SYS_config.debug_level>=1
        fprintf('              ');
    end
end
