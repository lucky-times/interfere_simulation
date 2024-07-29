function [eNodeB_sites,eNodeB_sectors,networkPathlossMap,networkShadowFadingMap] = NR_generate_network(SYS_config)
% �������˲����ļ�
% ���������
% SYS_config����������
% ���������
% eNodeB_sites��eNodeBʵ�弯
% eNodeB_sectors������ʵ�弯
% networkPathlossMap��·��ӳ�����
% networkShadowFadingMap����Ӱ˥��ӳ�����
%

if SYS_config.debug_level>=1 %���������Ϣ
    sprintf('Generating network\n');
end
%% �������˺����ó���
[eNodeB_sites,eNodeB_sectors,networkPathlossMap] = NR_init_network(SYS_config);

%% ���õ�Ƶ�ź��������źŵı�������ƽ̨�����ǵ�Ƶ�ź�
for b_=1:length(eNodeB_sites)
    for s_=1:length(eNodeB_sites(b_).sectors)
        data_power      = eNodeB_sites(b_).sectors(s_).max_power * (1-SYS_config.signaling_ratio);%signaling_ration=0
        signaling_power = eNodeB_sites(b_).sectors(s_).max_power * SYS_config.signaling_ratio;
        eNodeB_sites(b_).sectors(s_).max_power       = data_power;% �������ݹ��ʣ������ֻ�����
        eNodeB_sites(b_).sectors(s_).signaling_power = signaling_power;% ��Ƶ���ʣ�����Ϊ0
    end
end

%% ��ͬƵ����Ƶ��վ
num_hetnet_sites = networkPathlossMap.num_hetnet_sites;% �칹��վ��
num_first_sites = networkPathlossMap.num_first_sites;% ��һ��ϵͳ��վ��
for b_ = 1:length(eNodeB_sites)
    NR_init_get_eNodeB_neighbors(SYS_config,eNodeB_sites(b_), eNodeB_sites,num_hetnet_sites,num_first_sites);
end

%% ������Ӱ˥��
networkShadowFadingMap = NR_generate_shadowfading(networkPathlossMap,eNodeB_sites,SYS_config);

%% ��LTEϵͳ����MCL
% 5G�����᰸û��MCL
if strcmp(SYS_config.macroscopic_pathloss_model,'TS36942')
    % �����һ��ϵͳΪLTE
    if SYS_config.isDouble
        num_first_sectors = networkPathlossMap.num_first_sectors;
    else
        num_first_sectors = length(eNodeB_sectors);
    end
    LTE_sectors = 1:num_first_sectors;
    if strcmp(SYS_config.macroscopic_pathloss_model_settings.environment,'urban_macro')
        MCL = 70;% urban ����MCLΪ70dB
    elseif strcmp(SYS_config.macroscopic_pathloss_model_settings.environment,'rural_macro')
        MCL = 80;% rural ����MCLΪ80dB
    end
%     networkPathlossMap.apply_MCL(LTE_sectors,MCL);
end
if SYS_config.isDouble
    if strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')
        % ����ڶ���ϵͳΪLTE
        num_first_sectors = networkPathlossMap.num_first_sectors;
        LTE_sectors = num_first_sectors+1:length(eNodeB_sectors);
        if strcmp(SYS_config.macroscopic_pathloss_model_settings2.environment,'urban_macro')
            MCL = 70;
        elseif strcmp(SYS_config.macroscopic_pathloss_model_settings2.environment,'rural_macro')
            MCL = 80;
        end
        networkPathlossMap.apply_MCL(LTE_sectors,MCL);
    end
end

%% ������·����뷢�书�����վ���Ƿ�Χ
% NR_calculate_cell_coverage(SYS_config,networkPathlossMap,eNodeB_sites,eNodeB_sectors,networkShadowFadingMap);

% To avoid error of missing return argument
if ~exist('networkShadowFadingMap','var')
    networkShadowFadingMap = [];
end
