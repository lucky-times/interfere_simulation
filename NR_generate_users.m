function [ UEs ] = NR_generate_users(SYS_config,eNodeB_sites,eNodeB_sectors,networkPathlossMap,networkClock)
% ���㲢����UE
% ���������
% SYS_config����������
% eNodeB_sites����վ
% eNodeB_sectors������
% networkPathlossMap��·��ӳ�����
% networkClockʱ��ʵ��
%
% ���������
% UEs���û�
%

%% ����
% UE_spatial_distribution = spatial_distributions.constantUesPerCellUeSpatialDistribution(networkPathlossMap,eNodeB_sites,eNodeB_sectors,SYS_config);% ʵ����UE���㷽ʽ
% UE_positions = UE_spatial_distribution.generate_UE_positions(SYS_config);% UE_positions������UE��λ�þ������ǻ�û��UEʵ��

r = SYS_config.UE_r;
step = 0.5;
theta = 0:step:360-step;
UE_positions = [r*cosd(theta); r*sind(theta)]';

%% �����û�
UEs = network_elements.UE;
UE_pos_vector = zeros(size(UE_positions,1),2);
for u_ = 1:size(UE_positions,1)
    % �û���������
    UEs(u_)     = network_elements.UE; %ʵ�彨��
    UEs(u_).id  = u_; %�û�id����
%     UEs(u_).pos = NR_common_pixel_to_pos( UE_positions(u_,:), networkPathlossMap.coordinate_origin, networkPathlossMap.data_res); %�û�λ������
    UEs(u_).pos = UE_positions(u_,:);
    UE_pos_vector(u_,:) = UEs(u_).pos;
    UEs(u_).walking_model = walking_models.straightWalkingModel(SYS_config.UE_speed*SYS_config.TTI_length); % UE�˶�ģ��  %�뺯��ģ�岻ƥ�䣬ģ�����Ϊ���򣬴�ʱ����Ϊ·��
end
networkPathlossMap.UE_pos_vector = UE_pos_vector;

%% ��ʼ��RB����
% RB����ԭ����ƽ̨�̳й�����ֻ�õ���ʹ��RB��ȷ����վ���û�����Ĵ���
% ���Ҫ������Դ���䣬�����ʵ��޸����
NR_init_RB_grid(SYS_config,eNodeB_sites);

%% UE�������ֵĳ�ʼ��
for u_=1:length(UEs)
    
    % UE����
    antennas.antenna.attach_antenna_to_UE(UEs(u_),SYS_config);
    
    % ��������ϵ������
    UEs(u_).receiver_noise_figure = SYS_config.UE_receiver_noise_figure;
    UEs(u_).BS_receiver_noise_figure = SYS_config.BS_receiver_noise_figure;
    
    % UL power control error,�������������Ϊ����֪��ÿ��UE��power control error���
    UEs(u_).PCe = normrnd(0,SYS_config.PCe_param);
    %     UEs(u_).PCe = -9+18*rand(1,1);
    
    % ���չ�����������W��
    UEs(u_).thermal_noise_W_RB = 10^(0.1*SYS_config.UE.thermal_noise_density)/1000 * SYS_config.RB_bandwidth * 10^(UEs(u_).receiver_noise_figure/10);
    UEs(u_).BS_thermal_noise_W_RB = 10^(0.1*SYS_config.UE.thermal_noise_density)/1000 * SYS_config.RB_bandwidth * 10^(UEs(u_).BS_receiver_noise_figure/10);
    
    % ����ʱ��ģ��
    UEs(u_).clock = networkClock;
end

%% ȷ��С������UE
if length(UEs)==1 && isempty(UEs(1).id)
    no_UEs = true;
else
    no_UEs = false;
end
BS_pos  = zeros(length(eNodeB_sites), 2);
for i = 1:length(eNodeB_sites)
    BS_pos(i, :) = eNodeB_sites(i).pos;
end

if ~no_UEs
   % ����ǰ����õĻ�վ����������UE��λ�ã���UE��������վ
    UE_positions_m = zeros(length(UEs),2);
    
    num_hetnet_sectors = networkPathlossMap.num_hetnet_sectors;
    if SYS_config.isDouble
        switch SYS_config.scene_type
            case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
                num_first_UEs = (length(eNodeB_sectors)-num_hetnet_sectors)*SYS_config.UE_per_eNodeB;
            otherwise
                num_first_UEs = length(UEs)/2;
        end
    else
        num_first_UEs = 0;
    end
    networkPathlossMap.num_first_UEs = num_first_UEs;
    num_first_sectors = networkPathlossMap.num_first_sectors;
    for u_ = 1:length(UEs)     
        a = 1;
        % ���ݻ�վ����������������վ
        eNodeB_id = NR_calculate_attached_sector(SYS_config, BS_pos, UEs(u_).pos);
%         if theta(u_) <= 120
%            eNodeB_id = 31;
%         elseif theta(u_) <= 240
%            eNodeB_id = 32;
%         else
%            eNodeB_id = 33;
%         end
%         [ ~,~,eNodeB_id]    = networkPathlossMap.cell_assignment(SYS_config,u_,UEs(u_).pos,num_first_UEs);%��UE���ڵ�������eNB�ҳ���
        
       % ��UE�鸽����Ӧ��С��
        eNodeB_sectors(eNodeB_id).attachUser(UEs(u_));
        UE_positions_m(u_,:) = UEs(u_).pos;
        
        % ȷ��ͳ����ЩUE
        if SYS_config.isWraparound
            compute_only_UEs_from_this_eNodeBs = 1:57;% ֻͳ�Ʊ�1��19�Ż�վ��1��57��С���������UE����������
        else
            if SYS_config.isDouble
                if SYS_config.isS2F % �칹�д󳡾���С������С�������󳡾�
                    compute_only_UEs_from_this_eNodeBs = 1:num_first_sectors;% ��󳡾���С�����ĸ���
                else
                    compute_only_UEs_from_this_eNodeBs = num_first_sectors+1:length(eNodeB_sectors);% ��С�����Դ󳡾��ĸ���
                end
            else
                compute_only_UEs_from_this_eNodeBs = 1:length(eNodeB_sectors);
            end
        end
        if isempty(find(UEs(u_).attached_eNodeB.eNodeB_id == compute_only_UEs_from_this_eNodeBs,1))
            % ��ͳ�Ƹ�UE
            UEs(u_).deactivate_UE = true;%�ر�UEͳ��
        else
            % ͳ�Ƹ�UE
            UEs(u_).deactivate_UE = false;
        end
    end

end

function NR_init_RB_grid(SYS_config,eNodeB_sites)
% ��ʼ��RB

% Ϊÿ��С������RB
for b_ = 1:length(eNodeB_sites)
    for s_=1:length(eNodeB_sites(b_).sectors)
        
        % ����eNodeB�Ƿ�Ϊ����ģʽ�������Ƿ��з����û�
        eNodeB_sites(b_).sectors(s_).always_on = SYS_config.always_on;
        
        % RB �����Ĵ������ʼ��
        eNodeB_sites(b_).sectors(s_).RB_grid = network_elements.resourceBlockGrid(SYS_config.N_RB);
        eNodeB_sites(b_).sectors(s_).RB_grid.set_homogeneous_power_allocation(SYS_config,eNodeB_sites(b_).sectors(s_).max_power,eNodeB_sites(b_).sectors(s_).signaling_power);
    end
end
