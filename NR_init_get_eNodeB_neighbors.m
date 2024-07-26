function [in_interf_eNodeB_sites, in_interf_eNodeB_sectors, ad_interf_eNodeB_sites, ad_interf_eNodeB_sectors] = NR_init_get_eNodeB_neighbors(SYS_config,this_eNodeB,eNodeB_sites,num_hetnet_sites,num_first_sites)
% ��ĳ��վ��ͬƵ����Ƶ���Ż�վ
% �������: 
% SYS_config���������
% this_eNodeB����ǰ��վ
% eNodeB_sites�����л�վ
% num_hetnet_sites���칹�����У��ڶ���ϵͳ��վ��Ŀ
% num_first_sites����һ��ϵͳ��վ��
% �������:   
% in_interf_eNodeB_sites��ͬƵ���Ż�վ
% in_interf_eNodeB_sectors��ͬƵ��������
% ad_interf_eNodeB_sites����Ƶ���Ż�վ
% ad_interf_eNodeB_sectors����Ƶ��������

if SYS_config.isDouble
    switch  SYS_config.scene_type
        case {'UMa_to_UMi','UMa_to_InH','UMi_to_InH'}
            % ������UMa��wraparound��һ����75�У�
            if this_eNodeB.id>=1 && this_eNodeB.id<=(length(eNodeB_sites)-num_hetnet_sites)
                index = 1:(length(eNodeB_sites)-num_hetnet_sites);% ��һ��ϵͳȫ����վ
                index(find(index(:,:)==this_eNodeB.id)) = [];% �ų��Լ�
                in_interf_eNodeB_sites = eNodeB_sites(index);% ��ΪͬƵ��վ
                index = length(eNodeB_sites)-num_hetnet_sites+1:length(eNodeB_sites);% �ڶ���ϵͳ��ȫ����վ��Ϊ��Ƶ��վ
                ad_interf_eNodeB_sites = eNodeB_sites(index);
            else
                index = length(eNodeB_sites)-num_hetnet_sites+1:length(eNodeB_sites);
                index(find(index(:,:)==this_eNodeB.id)) = [];
                in_interf_eNodeB_sites = eNodeB_sites(index);
                index = 1:(length(eNodeB_sites)-num_hetnet_sites);
                ad_interf_eNodeB_sites = eNodeB_sites(index);
            end
            
        otherwise
            %% �һ�վ
            if length(eNodeB_sites)==1
                in_interf_eNodeB_sites = [];
                ad_interf_eNodeB_sites = [];
            else
                if SYS_config.isWraparound
                    if (this_eNodeB.id>=1 && this_eNodeB.id<=num_first_sites) % ��һ��ϵͳ
                        index_eNodeBs = 1:num_first_sites;% ��һ��ϵͳȫ����վ
                        index_eNodeBs(find(index_eNodeBs == this_eNodeB.id)) = []; % �ų����Լ�
                        allDistances = zeros(1,length(index_eNodeBs));
                        for b_ = 1:length(index_eNodeBs)
                            allDistances(b_) = sqrt((eNodeB_sites(index_eNodeBs(b_)).pos(1)-this_eNodeB.pos(1))^2+(eNodeB_sites(index_eNodeBs(b_)).pos(2)-this_eNodeB.pos(2))^2);
                        end
                        [~,index] = sort(allDistances);% ����index(1)=1ָindex_eNodeBs(1)���뱾��վ�����index(n)=5ָindex_eNodeBs(5)���뱾��վ�ŵ�n
                        lowest_idxs = index(1:18); % �ڵ�һ��ϵͳ���Ҿ����Լ������18����վ
                        lowest_idxs = sort(lowest_idxs);
                        index1 = index_eNodeBs(lowest_idxs);
                        in_interf_eNodeB_sites = eNodeB_sites(index1);
                        
                        index2 = [index1 this_eNodeB.id]+num_first_sites;
                        index2 = sort(index2);
                        ad_interf_eNodeB_sites = eNodeB_sites(index2);
                    else
                        % �ڶ���ϵͳ��վ
                        index_eNodeBs = num_first_sites+1:length(eNodeB_sites);
                        index_eNodeBs(find(index_eNodeBs == this_eNodeB.id)) = [];
                        allDistances = zeros(1,length(index_eNodeBs));
                        for b_ = 1:length(index_eNodeBs)
                            allDistances(b_) = sqrt((eNodeB_sites(index_eNodeBs(b_)).pos(1)-this_eNodeB.pos(1))^2+(eNodeB_sites(index_eNodeBs(b_)).pos(2)-this_eNodeB.pos(2))^2);
                        end
                        [~,index] = sort(allDistances);
                        lowest_idxs = index(1:18);
                        lowest_idxs = sort(lowest_idxs);
                        index1 = index_eNodeBs(lowest_idxs);
                        in_interf_eNodeB_sites = eNodeB_sites(index1);
                        
                        index2 = [index1 this_eNodeB.id]-num_first_sites;
                        index2 = sort(index2);
                        ad_interf_eNodeB_sites = eNodeB_sites(index2);
                    end
                else
                    if this_eNodeB.id>length(eNodeB_sites)/2 % �ڶ���ϵͳ�Ļ�վ
                        index = length(eNodeB_sites)/2+1:length(eNodeB_sites);% �ڶ���ϵͳ��ȫ����վ
                        index(find(index(:,:)==this_eNodeB.id)) = [];% Ȼ���ų��Լ�
                        in_interf_eNodeB_sites = eNodeB_sites(index);% ������ΪͬƵ��վ
                        index = 1:length(eNodeB_sites)/2;
                        ad_interf_eNodeB_sites = eNodeB_sites(index);% ��һϵͳ��ȫ����վ��Ϊ��Ƶ��վ
                    else % ��һ��ϵͳ�Ļ�վ
                        index = 1:length(eNodeB_sites)/2; % ��һ��ϵͳ��ȫ����վ
                        index(find(index(:,:)==this_eNodeB.id)) = []; % Ȼ���ų��Լ�
                        in_interf_eNodeB_sites = eNodeB_sites(index); % ������ΪͬƵ��վ
                        index = length(eNodeB_sites)/2+1:length(eNodeB_sites);
                        ad_interf_eNodeB_sites = eNodeB_sites(index); % �ڶ�ϵͳ��ȫ����վ��Ϊ��Ƶ��վ
                    end
                end
            end
    end
    
    this_eNodeB.in_interf_eNodeB_sites = in_interf_eNodeB_sites;
    this_eNodeB.ad_interf_eNodeB_sites = ad_interf_eNodeB_sites;
%             if this_eNodeB.id==62
%             fprintf('��ϵͳ���Ż�վ\n');
%             this_eNodeB.in_interf_eNodeB_sites.id
%             fprintf('*******\n');
%             fprintf('����ϵͳ���Ż�վ\n');
%             this_eNodeB.ad_interf_eNodeB_sites.id
%             fprintf('*******\n');
%             end
    
    %% ������
    if ~isempty(in_interf_eNodeB_sites)
        for s_=1:length(this_eNodeB.sectors)
            in_interf_eNodeB_sectors = this_eNodeB.sectors([1:(s_-1) (s_+1):end]); % ��ͬһ����վ�µ�ͬƵ��������
            
            for b_=1:length(in_interf_eNodeB_sites)
                in_interf_eNodeB_sectors = [in_interf_eNodeB_sectors in_interf_eNodeB_sites(b_).sectors]; % ��ͬ��վ�µ�ͬƵ��������
            end
            
            ad_interf_eNodeB_sectors=[];
            for b_=1:length(ad_interf_eNodeB_sites)
                ad_interf_eNodeB_sectors =[ad_interf_eNodeB_sectors ad_interf_eNodeB_sites(b_).sectors];
            end
            
            this_eNodeB.sectors(s_).in_interf_eNodeB_sectors = in_interf_eNodeB_sectors;
            this_eNodeB.sectors(s_).ad_interf_eNodeB_sectors = ad_interf_eNodeB_sectors;
            %                             this_eNodeB.id
            %                             this_eNodeB.sectors(s_).eNodeB_id
            %                             fprintf('��ϵͳ��������\n');
            %                             this_eNodeB.sectors(s_).neighbors_eNodeB.eNodeB_id
            %                             fprintf('*******\n');
            %                             fprintf('����ϵͳ��������\n');
            %                             this_eNodeB.sectors(s_).interf_eNodeB.eNodeB_id
            %                             fprintf('*******\n');
        end
    else
        for s_=1:length(this_eNodeB.sectors)
            this_eNodeB.sectors(s_).in_interf_eNodeB_sectors = [];
            this_eNodeB.sectors(s_).ad_interf_eNodeB_sectors = [];
        end
    end
    %һ��ϵͳ
else
    %% �һ�վ
    if length(eNodeB_sites)==1
        in_interf_eNodeB_sites = [];
    else
        if SYS_config.isWraparound
            allSites = reshape([eNodeB_sites.pos],2,[])';
            allDistances = sqrt(sum((repmat(this_eNodeB.pos,[size(allSites,1) 1])-allSites).^2,2));
            [allDistances_sorted indexes] = sortrows(allDistances);
            to_take = 19;
            lowest_idxs = indexes(1:min(length(indexes),to_take));% ȡ���Լ������19����վ�������Լ�
            lowest_idxs = sortrows(lowest_idxs);
            lowest_idxs(find(lowest_idxs == this_eNodeB.id)) = [];
            in_interf_eNodeB_sites = eNodeB_sites(lowest_idxs);
        else
            index = 1:length(eNodeB_sites);
            index(find(index(:,:)==this_eNodeB.id)) = [];
            in_interf_eNodeB_sites = eNodeB_sites(index);
        end
    end
    
    this_eNodeB.in_interf_eNodeB_sites = in_interf_eNodeB_sites;
    %         this_eNodeB.id
    %         fprintf('!!!!!!!!\n');
    %         this_eNodeB.neighbors_eNodeB.id
    %         fprintf('*******\n');
    
    %% ������
    if ~isempty(in_interf_eNodeB_sites)
        for s_=1:length(this_eNodeB.sectors)
            in_interf_eNodeB_sectors = this_eNodeB.sectors([1:(s_-1) (s_+1):end]); % ��ͬһ����վ�ڵ�ͬƵ��������
            
            for b_=1:length(in_interf_eNodeB_sites)
                in_interf_eNodeB_sectors = [in_interf_eNodeB_sectors in_interf_eNodeB_sites(b_).sectors];% ��ͬ��վ�µ�ͬƵ��������
            end
            
            this_eNodeB.sectors(s_).in_interf_eNodeB_sectors = in_interf_eNodeB_sectors;
        end
    else
        for s_=1:length(this_eNodeB.sectors)
            this_eNodeB.sectors(s_).in_interf_eNodeB_sectors = [];
        end
    end
end
end

