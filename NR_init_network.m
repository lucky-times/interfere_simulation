function [eNodeB_sites,eNodeB_sectors,networkMacroscopicPathlossMap] = NR_init_network(SYS_config)
% �������˲����ļ�
% ���������
% SYS_config����������ϵͳ
% ���������
% eNodeBs_sites����վ
% eNodeB_sectors������
% networkMacroscopicPathlossMap��·��ӳ��ͼ
%

%% �������ݸ�ֵ
data_res               = SYS_config.map_resolution;           % ��ͼ�ֱ���
eNodeB_sector_tx_power = SYS_config.eNodeB_tx_power;          % ��վ���书��

%% ����eNodeBʵ��
if SYS_config.debug_level>=1
    fprintf('Creating eNodeBs\n');
end
% ȷ����վ��λ��
[eNodeB_sites,num_hetnet_sites] = NR_init_create_eNodeBs(SYS_config);
% �ı�femto��վ��λ�ã����ݷ��������󷿼�����
if SYS_config.isFemto
    for b_ = 1:length(eNodeB_sites)
        if strcmp(eNodeB_sites(b_).site_type,'indoor')
            room_centre = eNodeB_sites(b_).pos;% ԭ����λ���Ƿ��������
            room_x = [room_centre(1)-10 room_centre(1)+10];
            room_y = [room_centre(2)-10 room_centre(2)+10];
            eNodeB_sites(b_).femto_room = [room_x;room_y];% ���淿�����򣬻�ͼ��ʱ����
            new_pos = [random('unif',room_x(1),room_x(2)),random('unif',room_y(1),room_y(2))];% ��femto��վһ�������λ��
            eNodeB_sites(b_).pos = new_pos;
        end
    end
end

% ����վ����������ֵ
s_idx   = 1;
eNodeB_sectors = network_elements.eNodeB_sector;  % ��ʼ��,eNodeBsΪeNodeB_sector�Ķ���

for b_ = 1:length(eNodeB_sites)-num_hetnet_sites
   % ����С��ʵ��
    eNodeB_sites(b_).sectors    = network_elements.eNodeB_sector;
    for s_ = 1:length(SYS_config.sector_azimuths)
        eNodeB_sites(b_).sectors(s_)               = network_elements.eNodeB_sector;
        eNodeB_sites(b_).sectors(s_).parent_eNodeB = eNodeB_sites(b_);
        eNodeB_sites(b_).sectors(s_).id            = s_; %����С��id������ڱ�eNodeB���ԣ�ȡֵΪ1,2,3��
        switch SYS_config.scene_type
            case {'UMA','RMa'}
                eNodeB_sites(b_).sectors(s_).azimuth       = utils.miscUtils.wrapTo359(SYS_config.antenna_azimuth_offsett + SYS_config.sector_azimuths(s_));
            case {'UMI','UMa_to_UMi'}
                if ~SYS_config.isManhattan %Rand drop ������ָ�򸲸�ȦԲ��
                    eNodeB_sites(b_).sectors(s_).azimuth = atan2(eNodeB_sites(b_).parent_centre_pos(2) -eNodeB_sites(b_).pos(2),eNodeB_sites(b_).parent_centre_pos(1)-eNodeB_sites(b_).pos(1))./pi*180;
                else
                    eNodeB_sites(b_).sectors(s_).azimuth = utils.miscUtils.wrapTo359(SYS_config.antenna_azimuth_offsett + 330);
                end
            case {'InH','UMa_to_InH','UMi_to_InH','InH2','InH3'} %InH������ƫ��Ϊ0
                eNodeB_sites(b_).sectors(s_).azimuth = 0;
        end
        eNodeB_sites(b_).sectors(s_).max_power     = eNodeB_sector_tx_power; %���û�վ������书��
        eNodeB_sites(b_).sectors(s_).antenna_type  = SYS_config.antenna.antenna_gain_pattern; % ��������
        
        eNodeB_sites(b_).sectors(s_).eNodeB_id     = s_idx; %����С����ȫ��id��
        eNodeB_sectors(s_idx) = eNodeB_sites(b_).sectors(s_); %��������ͨ��eNodeB_sites(b_).sectors(s_)�������ģ����ڽ��丳��eNodeB_sectors���Ժ����ͨ��eNodeB_sectors���ж������Ĳ���
        
         % Ϊ��վ��������
        antennas.antenna.attach_antenna_to_eNodeB(eNodeB_sites(b_).sectors(s_),SYS_config,eNodeB_sites(b_).site_type);% ��ʼ������
        
        % ������·��ģ��
        if SYS_config.debug_level>=1
            fprintf('Site %d, eNodeB %d: ',b_,s_);%b_Ϊ��վ��s_Ϊ����
        end
        
        %˫ϵͳ
        if SYS_config.isDouble 
            switch SYS_config.scene_type
                case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
                    % ���칹��˵���ⲽ���ǵ�һϵͳ��
                    eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
                        SYS_config,...
                        SYS_config.macroscopic_pathloss_model,...
                        SYS_config.frequency,...
                        SYS_config.macroscopic_pathloss_model_settings);
                otherwise
                    if b_<=length(eNodeB_sites)/2 %��һ��ϵͳ·��ģ��
                        eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
                            SYS_config,...
                            SYS_config.macroscopic_pathloss_model,...
                            SYS_config.frequency,...
                            SYS_config.macroscopic_pathloss_model_settings);
                    else %�ڶ���ϵͳ·��ģ��
                        eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
                            SYS_config,...
                            SYS_config.macroscopic_pathloss_model2,...
                            SYS_config.frequency2,...
                            SYS_config.macroscopic_pathloss_model_settings2);
                    end
            end
        else %��ϵͳ��·��ģ��
            eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
                SYS_config,...
                SYS_config.macroscopic_pathloss_model,...
                SYS_config.frequency,...
                SYS_config.macroscopic_pathloss_model_settings);
        end
        s_idx = s_idx + 1;
    end
end
num_hetnet_sectors = 0;% �칹�����У��ڶ���ϵͳ��������
for b_ = length(eNodeB_sites)-num_hetnet_sites+1:length(eNodeB_sites)% �������칹������β�ִ��,num_hetnet_sites~=0�Ļ��������칹����
    eNodeB_sites(b_).sectors = network_elements.eNodeB_sector;
    for s_ = 1:length(SYS_config.sector_azimuths2)
        eNodeB_sites(b_).sectors(s_)               = network_elements.eNodeB_sector;
        eNodeB_sites(b_).sectors(s_).parent_eNodeB = eNodeB_sites(b_);
        eNodeB_sites(b_).sectors(s_).id            = s_;      
        switch eNodeB_sites(b_).site_type
            case 'macro'
                eNodeB_sites(b_).sectors(s_).azimuth = utils.miscUtils.wrapTo359(SYS_config.antenna_azimuth_offsett + SYS_config.sector_azimuths2(s_));
            case 'micro'
                eNodeB_sites(b_).sectors(s_).azimuth = atan2(eNodeB_sites(b_).parent_centre_pos(2) -eNodeB_sites(b_).pos(2),eNodeB_sites(b_).parent_centre_pos(1)-eNodeB_sites(b_).pos(1))./pi*180;
        end
        eNodeB_sites(b_).sectors(s_).max_power     = SYS_config.eNodeB_tx_power2;
        eNodeB_sites(b_).sectors(s_).antenna_type  = SYS_config.antenna.antenna_gain_pattern;
        eNodeB_sites(b_).sectors(s_).eNodeB_id     = s_idx;
        eNodeB_sectors(s_idx)                              = eNodeB_sites(b_).sectors(s_);
        antennas.antenna.attach_antenna_to_eNodeB(eNodeB_sites(b_).sectors(s_),SYS_config,eNodeB_sites(b_).site_type);
        
        if SYS_config.debug_level>=1
            fprintf('Site %d, eNodeB %d: ',b_,s_);%b_Ϊ��վ��s_Ϊ����
        end
        
        eNodeB_sites(b_).sectors(s_).macroscopic_pathloss_model = macroscopic_pathloss_models.generalPathlossModel.generateMacroscopicPathlossModel(...
            SYS_config,...
            SYS_config.macroscopic_pathloss_model2,...
            SYS_config.frequency2,...
            SYS_config.macroscopic_pathloss_model_settings2);
        
        s_idx = s_idx + 1;
        num_hetnet_sectors = num_hetnet_sectors + 1;
    end
end

%% ȷ����ͼ��С
switch SYS_config.scene_type
    case {'UMA','RMa','UMa_to_InH'}
        % �õ���վ��λ��
        tx_pos = zeros(length(eNodeB_sites),2);%2�У���һ��Ϊx���ڶ���Ϊy
        for b_ = 1:length(eNodeB_sites)
            tx_pos(b_,:) = eNodeB_sites(b_).pos;
        end
        
        % ���ݻ�վ��λ�ü����ͼ�߽�
        if  SYS_config.nr_eNodeB_rings ~= 0
            roi_x = [min(tx_pos(:,1)),max(tx_pos(:,1))];
            roi_y = [min(tx_pos(:,2)),max(tx_pos(:,2))];
        else % �����վֻ��һ��
            roi_x = [-SYS_config.ISD,SYS_config.ISD];
            roi_y = [-SYS_config.ISD,SYS_config.ISD];
        end
        
        ISD = SYS_config.ISD;
        roi_x =[roi_x(1)-2*ISD/3 roi_x(2)+2*ISD/3];
        roi_y =[roi_y(1)-2*ISD/3 roi_y(2)+2*ISD/3];
    case 'UMI'
        if ~SYS_config.isManhattan
            % ����ISD������ģ��ܸպð���ȫ����վ����վ�ĸ��Ƿ�Χ
            ISD = SYS_config.ISD;
            roi_x =[-(ISD*2)-(ISD*2)/3 (ISD*2)+(ISD*2)/3];
            roi_y =[-(ISD*(3^0.5))-(ISD*2)/3 (ISD*(3^0.5))+(ISD*2)/3];
        else
            % �̶��Ĵ�С
            roi_x =[0 1080];
            roi_y =[0 1080];
        end
    case 'UMa_to_UMi'
        if ~SYS_config.isManhattan
            ISD = SYS_config.ISD;
            roi_x =[-(ISD*2)-(ISD*2)/3 (ISD*2)+(ISD*2)/3];
            roi_y =[-(ISD*(3^0.5))-(ISD*2)/3 (ISD*(3^0.5))+(ISD*2)/3];
        else
            tx_pos = zeros(length(eNodeB_sites),2);
            for b_ = 1:length(eNodeB_sites)
                tx_pos(b_,:) = eNodeB_sites(b_).pos;
            end
            roi_x = [min(tx_pos(:,1)),max(tx_pos(:,1))];
            roi_y = [min(tx_pos(:,2)),max(tx_pos(:,2))];
            ISD = SYS_config.ISD;
            roi_x =[roi_x(1)-2*ISD/3 roi_x(2)+2*ISD/3];
            roi_y =[roi_y(1)-2*ISD/3 roi_y(2)+2*ISD/3];
        end
    case 'UMi_to_InH'
        UMi_r = SYS_config.UMi_r;
        tx_pos = zeros(length(eNodeB_sites),2);
        for b_ = 1:length(eNodeB_sites)
            tx_pos(b_,:) = eNodeB_sites(b_).pos;
        end
        roi_x = [min(tx_pos(:,1)),max(tx_pos(:,1))];
        roi_y = [min(tx_pos(:,2)),max(tx_pos(:,2))];
        roi_x =[roi_x(1)-UMi_r roi_x(2)+UMi_r];
        roi_y =[roi_y(1)-UMi_r roi_y(2)+UMi_r];
    case {'InH','InH2'}
        % �̶��Ĵ�С
        roi_x =[0 120];
        if ~SYS_config.isFemto
            roi_y =[0 50];
        else
            roi_y =[0 40];
        end
    case {'InH3'}
        % �̶��Ĵ�С
        roi_x =[0 120];
        if ~SYS_config.isFemto
            roi_y =[0 100+SYS_config.InH3_d];
        else
            roi_y =[0 40];
        end
end

%% ��ʼ��·��ͼ��
networkMacroscopicPathlossMap                        = channel_gain_wrappers.macroscopicPathlossMap;%�޲ι��캯��
networkMacroscopicPathlossMap.data_res               = data_res;
networkMacroscopicPathlossMap.roi_x                  = roi_x;
networkMacroscopicPathlossMap.roi_y                  = roi_y;
networkMacroscopicPathlossMap.num_hetnet_sites       = num_hetnet_sites;
networkMacroscopicPathlossMap.num_hetnet_sectors     = num_hetnet_sectors;

if SYS_config.isDouble
    switch SYS_config.scene_type
        case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
            num_first_sites = length(eNodeB_sites) - num_hetnet_sites;
            num_first_sectors = length(eNodeB_sectors)-num_hetnet_sectors;
        otherwise
            num_first_sites = length(eNodeB_sites)/2;
            num_first_sectors = length(eNodeB_sectors)/2;
    end
else
    % ��ϵͳͨ������length��eNodeB_sites�����ýϿ�
    num_first_sites = 0;
    num_first_sectors = 0;
end
networkMacroscopicPathlossMap.num_first_sectors = num_first_sectors;% ˫ϵͳ�е�һ��ϵͳ��������
networkMacroscopicPathlossMap.num_first_sites = num_first_sites;% ˫ϵͳ�е�һ��ϵͳ�Ļ�վ��


%% ����·�����
if SYS_config.debug_level>=1
    fprintf('Creating cell pathloss map\n');
end
% macroscopic_pathloss_models.generalPathlossModel.calculate_pathloss_maps(SYS_config,eNodeB_sites,networkMacroscopicPathlossMap);

%% ����ǰ��ƽ̨�̳����ģ�ûʲô��
networkMacroscopicPathlossMap.name = 'NR or LTE pathloss model';
