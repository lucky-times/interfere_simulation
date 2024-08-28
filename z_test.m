%test.m
%���ļ������޸Ŀɱ����������ƽ̨����ͳ�����ݡ�

%��֮ǰ���е�ͼ���Լ�����ɾ���������ڹ۲������
close all force;
clc;
clear all
clear global;
clear classes;

SYS_config = NR_load_params;%����Ĭ�Ϸ������

%% ��Ĭ�ϲ����Ļ����϶�һЩ�ɱ����Ȥ���������޸ĵ���
% һЩ���õĿɱ�������£�
SYS_config.shadow_fading_type         = 'claussen';
SYS_config.compact_results_file       = true;
SYS_config.simulation_time_tti        = 1; %���澭����ʱ���
SYS_config.keep_UEs_still             = true; %UE�Ƿ񱣳־�ֹ
SYS_config.UE_speed = 5/3.6;% UE�˶��ٶ�5km/h,(5/3.6)m/s
SYS_config.TTI_length = 1;% Length of a TTI (subframe), in seconds.

% ������׼�ľ�ȷ��
SYS_config.beam_accuracy_phi_ue = 0; %UE��������ˮƽ�ǣ��ȣ�
SYS_config.beam_accuracy_theta_ue = 0; %UE�������ʹ�ֱ�ǣ��ȣ�
SYS_config.beam_accuracy_phi_bs = 0; %BS��������ˮƽ�ǣ��ȣ�
SYS_config.beam_accuracy_theta_bs = 0; %BS�������ʹ�ֱ�ǣ��ȣ�
SYS_config.beam_accurancy_type = 'Constant'; %Constant��ʾ��׼ƫ��Ϊ����ֵ������Unique��ʾ��׼ƫ��Ϊ0-����ֵ�ľ��ȷֲ�
SYS_config.asynchronization = false; %�����Ƿ�ͬ��

SYS_config.default_shown_GUI_cells = [];
SYS_config.frequency = 4.9e9; %victimϵͳ��Ƶ��

SYS_config.photo = true;
SYS_config.PCe_param = 0; %���صľ�ȷ��

SYS_config.isDouble = false;% �Ƿ�˫ϵͳ
SYS_config.interference_type = 0;% 0ͬƵ���ţ�1��Ƶ���ţ�2���У�ֻ��˫ϵͳ����Ч
SYS_config.antenna_mode=0;  % 0��ʾBS��UE���޲������ͣ�1��ʾBS�������ͣ�UE�޲������ͣ�2��ʾBS��UE����������
SYS_config.attatch_mode=3;  % 0��ʾ������·��1��ʾ·���BS�˵�Ԫ���棻2��ʾ·���UE�˵�Ԫ���棻3��ʾ·���BS��UE�˵�Ԫ���档����ȷ��UE����С�������㣩
SYS_config.use_cache = false;% ���������Ƿ�ʹ���Ѿ��洢�õ�����
SYS_config.isSave = false; %�Ƿ񱣴沨�����Ͳ���������
SYS_config.antenna_element = [8,16];% UMi UMaԭ����8,16
SYS_config.antenna_element_InH = [2,2];% InH 30Gԭ����4,8
SYS_config.antenna_theta_3db = 6;%��

SYS_config.bandwidth = 100e6; %����ϵͳ�Ĵ���
SYS_config.UE_tx_power = 10^2.3/1000; %UE�ķ��书�ʣ���λW���˴�Ϊ23dBm��
SYS_config.Gama = 1; %�����е�gamma����
SYS_config.n_UE_served_per_BS = 1; % ÿ��BSͬһʱ������UE��

SYS_config.antenna.antenna_gain_pattern = 'NRAntennaBeamforming'; %�������͵����ͣ������ڴ˴�������չ�����Ĳ�������ͼ����
SYS_config.UE_height = 200;%UAV���и߶� 100 200 300 500
SYS_config.UE_r = 800; %UAV�����Ļ�վ�ļ���  100 300 500 800 1200
str = sprintf("h=%d, r=%d", SYS_config.UE_height, SYS_config.UE_r);
disp(str)
% SYS_config.eNodeB_pos = zeros(19,2);
SYS_config.UE_max_antenna_gain = 5; %UE������������棨��Ԫ���棩
SYS_config.AntennaPattern3d = true;
SYS_config.UE_receiver_noise_figure = 10; %UE�˽��ջ�������ϵ��
SYS_config.BS_receiver_noise_figure = 10; %BS�˽��ջ�������ϵ��
SYS_config.macroscopic_pathloss_model = 'TS36942'; %��һ��ϵͳ·��ģ�͵�ѡ��TS38900��ʾNRϵͳ��TS36942��ʾLTE·��
SYS_config.macroscopic_pathloss_model2 = 'TS36942'; %�ڶ���ϵͳ·��ģ�͵�ѡ��TS38900��ʾNRϵͳ��TS36942��ʾLTE·��

SYS_config.asynchronization_switch=0;%�жϸ������ͣ�0��ʾͬ����1��ʾ�첽��2��ʾ��50%

SYS_config.cable_loss=0; %�����������
SYS_config.scene_type = 'UMA'; %����������ѡ��
SYS_config.beam_loss=0; %�������
SYS_config.tilt = 0; %��վ���������

switch SYS_config.scene_type
    case 'UMA' %UMA������һЩ�������趨
        SYS_config.isWraparound = false; %�Ƿ�ʹ��Wraparound
        SYS_config.isNewUMA = true;% ��ϵ���������û�����
        if SYS_config.isNewUMA % ѡ��ISD���°汾200 ISD���������û�����Ϊ0��10���ϰ汾500ISD���������û�����Ϊ8:2
            SYS_config.ISD = 1000; 
        else
            SYS_config.ISD = 500;
        end
        SYS_config.shift_mode = 0;%�Ƿ�ַ,0����ַ��1����ƫ��100%,2�Լ�����ƫ����
        if SYS_config.shift_mode == 2
            SYS_config.angle_between_systems = 45; %ƫ�ýǶȣ���λ��
            SYS_config.ISD_between_systems = (2^0.5)/2*SYS_config.ISD; %ƫ�ü������λ��
        end
        SYS_config.hand_over = false; %�Ƿ����3dBhandover
        SYS_config.sector_azimuths = 0:360/3:359; %�������ߵķ����
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_macro'; %��һ��ϵͳ�Ļ���
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro'; %�ڶ���ϵͳ�Ļ���
        if SYS_config.ISD <= 300
            if  SYS_config.isWraparound
                SYS_config.map_resolution = 10; %���ú�·��ͼ���ص������
                SYS_config.shadow_fading_map_resolution = 10;%������Ӱ˥��ͼ�����ص������
            else
                SYS_config.map_resolution = 5;
                SYS_config.shadow_fading_map_resolution=5;
            end
        else
            if  SYS_config.isWraparound
                SYS_config.map_resolution = 20;
                SYS_config.shadow_fading_map_resolution = 20;
            else
                SYS_config.map_resolution = 50;
                SYS_config.shadow_fading_map_resolution=10;
            end
        end
        SYS_config.n_snapshot = 1;%���ÿ�����
        SYS_config.shadow_fading_sd_LOS = 4; %LOS����Ӱ˥���׼��
        SYS_config.shadow_fading_sd_NLOS = 6; %NLOS����Ӱ˥���׼��
        SYS_config.r_eNodeBs = .5; %��Ӱ˥��վ�����ϵ��
        SYS_config.site_height=25; %��վ��վַ�߶�
        SYS_config.eNodeB_tx_power = 10^4.3/1000; %��վ������书�� w
        SYS_config.antenna.max_antenna_gain = 8; %BS�������������
        SYS_config.default_shown_GUI_cells = 1:57;
    case 'UMI' %����Random drop ��manhattan�����ӳ���
        SYS_config.ISD = 200;
        if SYS_config.ISD <= 300
            SYS_config.map_resolution = 5;
            SYS_config.shadow_fading_map_resolution=5;
        else
            SYS_config.map_resolution = 10;
            SYS_config.shadow_fading_map_resolution=10;
        end
        SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
        SYS_config.shift_mode = 0;%�Ƿ�ַ,1������ַ
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_micro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_micro';
        SYS_config.n_snapshot = 5;%���ÿ�����
        SYS_config.shadow_fading_sd_LOS = 4;
        SYS_config.shadow_fading_sd_NLOS = 7.82;
        SYS_config.r_eNodeBs = .5;
        SYS_config.site_height=10;
        SYS_config.eNodeB_tx_power = 10^3/1000;% ��λw
        SYS_config.isManhattan = false; %�Ƿ�������������ˣ�true�Ļ��������������ˣ��������Random drop����
        if SYS_config.isManhattan
            SYS_config.map_resolution = 5;
            SYS_config.shadow_fading_map_resolution=5;
        end
        SYS_config.antenna.max_antenna_gain = 8;
        SYS_config.default_shown_GUI_cells = [];
    case 'InH' %����Open office������Home femto����
        SYS_config.shift_mode = 0;%�Ƿ�ַ,0����ַ
        SYS_config.isFemto = false; %�Ƿ�Ϊfemto��������Ϊtrue��ΪHome femto����������ΪOpen office����
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
        SYS_config.map_resolution = 0.5;
        SYS_config.shadow_fading_map_resolution=0.5;
        SYS_config.n_snapshot = 10;%���ÿ�����
        SYS_config.shadow_fading_sd_LOS = 3;
        SYS_config.shadow_fading_sd_NLOS = 8.03;
        SYS_config.r_eNodeBs = 0;
        SYS_config.site_height=3;
        SYS_config.UE_height=1;%UE�ĸ߶�Ϊ1m
        SYS_config.eNodeB_tx_power = 10^2.3/1000;% ��λw     
        SYS_config.antenna.max_antenna_gain = 5;
        SYS_config.default_shown_GUI_cells = [];
    case 'RMa'
        SYS_config.isWraparound = true; %�Ƿ�ʹ��Wraparound
        SYS_config.frequency = 7e9;
        SYS_config.frequency2 = 7e9;
        SYS_config.ISD = 1732;
        SYS_config.shift_mode = 0;%�Ƿ�ַ,0����ַ��1����ƫ��100%,2�Լ�����ƫ����
        if SYS_config.shift_mode == 2
            SYS_config.angle_between_systems = 45;
            SYS_config.ISD_between_systems = 50;
        end
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0:360/3:359;
        SYS_config.macroscopic_pathloss_model_settings.environment = 'rural_macro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'rural_macro';
        if SYS_config.isWraparound
            SYS_config.map_resolution = 80;
            SYS_config.shadow_fading_map_resolution = 80;
        else
            SYS_config.map_resolution = 50;
            SYS_config.shadow_fading_map_resolution = 50;
        end
        SYS_config.n_snapshot = 10;%���ÿ�����
        SYS_config.shadow_fading_sd_LOS = 4;
        SYS_config.shadow_fading_sd_LOS2 = 6; %900��RMA������LOS״̬������LOS��
        SYS_config.shadow_fading_sd_NLOS = 8;
        SYS_config.r_eNodeBs = .5;
        SYS_config.site_height=35;
        SYS_config.eNodeB_tx_power = 10^4.3/1000;% ��λw
        SYS_config.antenna.max_antenna_gain = 8;
        SYS_config.default_shown_GUI_cells = 1:57;
    case 'UMa_to_UMi'
        SYS_config.isDouble = true;% �Ƿ�˫ϵͳ
        SYS_config.isNewUMA = true;
        SYS_config.isManhattan = false; %΢С���Ƿ��������٣�֧��UMa��UMi�����ӳ����Ĺ����о�
        SYS_config.isS2F = true;
        SYS_config.ISD = 200;
        if SYS_config.ISD <= 300
            SYS_config.map_resolution = 5;
            SYS_config.shadow_fading_map_resolution=5;
        else
            SYS_config.map_resolution = 10;
            SYS_config.shadow_fading_map_resolution=10;
        end
        SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
        SYS_config.shift_mode = 1;%�Ƿ�ַ,1������ַ
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;% ȫ������
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_micro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro';
        SYS_config.n_snapshot = 5;%���ÿ�����
        SYS_config.n_snapshot2 = 5;%���õڶ���ϵͳ������
        SYS_config.shadow_fading_sd_LOS = 4;
        SYS_config.shadow_fading_sd_NLOS = 7.82;
        SYS_config.r_eNodeBs = .5;
        SYS_config.site_height=10;%10
        SYS_config.eNodeB_tx_power = 10^3.3/1000;% ��λw
        SYS_config.antenna.max_antenna_gain = 8;
        SYS_config.default_shown_GUI_cells = [];
        
        if SYS_config.isManhattan
            SYS_config.shift_mode = 1;% ����ϵͳ���ƫ�ƣ�1����Ĭ�ϣ�������λ���м䣩��2�������ã��ڶ���ϵͳ��������������ĵ�λ�ã�
            if SYS_config.shift_mode == 2
                SYS_config.angle_between_systems = 45; %ƫ�ýǶ�
                SYS_config.ISD_between_systems = 200; %ƫ�ü��
            end
        end
        
        SYS_config.sector_azimuths2 = 0:360/3:359; %�칹����ϵͳ�����������
        SYS_config.eNodeB_tx_power2 = 10^4.3/1000; %�칹����ϵͳ��վ�ķ��书��
        SYS_config.site_height2 = 25;              %�칹����ϵͳ��վ�ĸ߶�
        SYS_config.shadow_fading_sd_LOS_het = 4;   %�칹����ϵͳ��վLOS�µ���Ӱ˥���׼��
        SYS_config.shadow_fading_sd_NLOS_het = 6;  %�칹����ϵͳ��վNLOS�µ���Ӱ˥���׼��
        SYS_config.antenna.max_antenna_gain2 = 8;  %�칹����ϵͳ��վ�������Ԫ����
    case 'UMa_to_InH'
        SYS_config.isDouble = true;% �Ƿ�˫ϵͳ
        SYS_config.isFemto = true;% ֧��UMa��InH�����ӳ������칹�����о�
        SYS_config.isS2F = true;
        SYS_config.ISD = 200;
        if SYS_config.ISD <= 300
            SYS_config.map_resolution = 2;
            SYS_config.shadow_fading_map_resolution=2;
        else
            SYS_config.map_resolution = 5;
            SYS_config.shadow_fading_map_resolution=5;
        end
        SYS_config.shift_mode = 1;% ����ϵͳ���ƫ�ƣ�1����Ĭ�ϣ�inhλ���м䣩��2�������ã��ڶ���ϵͳ�����inh���ĵ�λ�ã�
        if SYS_config.shift_mode == 2
            SYS_config.angle_between_systems = 45; %ƫ�ýǶ�
            SYS_config.ISD_between_systems = 30; %ƫ�ü��
        end
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;% ȫ������
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro';
        SYS_config.n_snapshot = 5;%���ÿ�����
        SYS_config.n_snapshot2 = 5;%���õڶ���ϵͳ������
        SYS_config.shadow_fading_sd_LOS = 3;
        SYS_config.shadow_fading_sd_NLOS = 8.03;
        SYS_config.r_eNodeBs = .5;
        SYS_config.site_height=3;
        SYS_config.UE_height=1.5;
        SYS_config.eNodeB_tx_power = 10^2.3/1000;% ��λw
        SYS_config.antenna.max_antenna_gain = 5;
        SYS_config.default_shown_GUI_cells = [];
        
        SYS_config.sector_azimuths2 = 0:360/3:359;
        SYS_config.eNodeB_tx_power2 = 10^4.3/1000;
        SYS_config.site_height2 = 25;
        SYS_config.shadow_fading_sd_LOS_het = 4;
        SYS_config.shadow_fading_sd_NLOS_het = 6;
        SYS_config.antenna.max_antenna_gain2 = 8;
    case 'UMi_to_InH'
        SYS_config.isDouble = true;% �Ƿ�˫ϵͳ
        SYS_config.isFemto = false;% ��֧��Random drop��InH�����ӳ������칹�����о�
        SYS_config.isS2F = true;
        SYS_config.ISD = 200;
        if SYS_config.ISD <= 300
            SYS_config.map_resolution = 2;
            SYS_config.shadow_fading_map_resolution=2;
        else
            SYS_config.map_resolution = 5;
            SYS_config.shadow_fading_map_resolution=5;
        end
        SYS_config.shift_mode = 1;% ����ϵͳ���ƫ�ƣ�1����Ĭ�ϣ�inhλ���м䣩��2�������ã��ڶ���ϵͳ�����inh���ĵ�λ�ã�
        if SYS_config.shift_mode == 2
            SYS_config.angle_between_systems = 45; %ƫ�ýǶ�
            SYS_config.ISD_between_systems = 30; %ƫ�ü��
        end
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;% ȫ������
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_micro';
        SYS_config.n_snapshot = 5;%���ÿ�����
        SYS_config.n_snapshot2 = 5;%���õڶ���ϵͳ������
        SYS_config.shadow_fading_sd_LOS = 3;
        SYS_config.shadow_fading_sd_NLOS = 8.03;
        SYS_config.r_eNodeBs = .5;
        SYS_config.site_height=3;
        SYS_config.UE_height=1.5;
        SYS_config.eNodeB_tx_power = 10^2.3/1000;% ��λw
        SYS_config.antenna.max_antenna_gain = 5;
        SYS_config.default_shown_GUI_cells = [];
        
        SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
        SYS_config.sector_azimuths2 = 0;
        SYS_config.eNodeB_tx_power2 = 10^3.3/1000;
        SYS_config.site_height2 = 10;
        SYS_config.shadow_fading_sd_LOS_het = 4;
        SYS_config.shadow_fading_sd_NLOS_het = 7.82;
        SYS_config.antenna.max_antenna_gain2 = 8;
        
    case 'InH2' %����Open office������Home femto����
        SYS_config.shift_mode = 0;%�Ƿ�ַ,0����ַ
        SYS_config.isFemto = false; %�Ƿ�Ϊfemto��������Ϊtrue��ΪHome femto����������ΪOpen office����
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
        SYS_config.map_resolution = 0.5;
        SYS_config.shadow_fading_map_resolution=0.5;
        SYS_config.n_snapshot = 10;%���ÿ�����
        SYS_config.shadow_fading_sd_LOS = 3;
        SYS_config.shadow_fading_sd_NLOS = 8.03;
        SYS_config.r_eNodeBs = 0;
        SYS_config.site_height=3;
        SYS_config.UE_height=1;%UE�ĸ߶�Ϊ1m
        SYS_config.eNodeB_tx_power = 10^2.3/1000;% ��λw     
        SYS_config.antenna.max_antenna_gain = 5;
        SYS_config.default_shown_GUI_cells = [];
        
    case 'InH3' %����Open office������Home femto����
        SYS_config.InH3_d = 10;
        SYS_config.shift_mode = 0;%�Ƿ�ַ,0����ַ
        SYS_config.isFemto = false; %�Ƿ�Ϊfemto��������Ϊtrue��ΪHome femto����������ΪOpen office����
        SYS_config.hand_over = true;
        SYS_config.sector_azimuths = 0;
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
        SYS_config.map_resolution = 0.5;
        SYS_config.shadow_fading_map_resolution=0.5;
        SYS_config.n_snapshot = 50;%���ÿ�����
        SYS_config.shadow_fading_sd_LOS = 3;
        SYS_config.shadow_fading_sd_NLOS = 8.03;
        SYS_config.r_eNodeBs = 0;macroscopic_pathloss_model_settings
        SYS_config.site_height=3;
        SYS_config.UE_height=1;%UE�ĸ߶�Ϊ1m
        SYS_config.eNodeB_tx_power = 10^2.3/1000;% ��λw     
        SYS_config.antenna.max_antenna_gain = 5;
        SYS_config.default_shown_GUI_cells = [];
end

output_results_file = NR_sim_main(SYS_config); %����������
%simulation_data                   = load(output_results_file); %��������
close all;

% GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data); %��������ͼ