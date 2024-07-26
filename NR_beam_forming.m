function [beam_antenna_gain,beam_ue_gain] = NR_beam_forming(UEs,eNodeBs_sectors,networkPathlossMap,SYS_config)
% ������������ͼ��
% ���������
% UEs���û��б�
% eNodeBs_sectors��С��ʵ��
% networkPathlossMap��·��ӳ�����
% SYS_config��LTE����ϵͳ
%
% ���������
% beam_antenna_gain��BS�˲�����������
% beam_ue_gain��UE�˲�����������

num_markers = 10; %ͬ�������������������
u_markings  = round(linspace(1,length(UEs),num_markers));
DEBUG_LEVEL = SYS_config.debug_level;

ue_len=length(UEs); %���е�UE����
enb_len=length(eNodeBs_sectors);%����ϵͳ�Ļ������122����վ������UMA wrap around��,������sector����ֵΪ366
beam_antenna_gain=zeros(ue_len,enb_len );
beam_ue_gain=zeros(ue_len,enb_len );
flag=true; %�����ж�ÿ��sector���û����Ƿ���ͬ
if DEBUG_LEVEL>=1
    if SYS_config.antenna_mode==2
        fprintf('Creating BS and UE beam forming: ');
    elseif SYS_config.antenna_mode==1
        fprintf('Creating BS beam forming: ');
    end
end


if SYS_config.antenna_mode==0 %ȷ��BS����UE���Ƿ�����˲�������
    if DEBUG_LEVEL>=1
        fprintf('No BF in BS nor UE\n');
    end
else
    %���벨����������
    if SYS_config.use_cache
        try
            filename1 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain.mat'];
            filename2 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain.mat'];
            tmp = load(filename1);
            beam_antenna_gain = tmp.beam_antenna_gain;
            tmp = load(filename2);
            beam_ue_gain = tmp.beam_ue_gain;
            creat_beam = false;
            if DEBUG_LEVEL>=1
                fprintf('(Use beam forming cache)\n');
            end
        catch
            if DEBUG_LEVEL>=1
                fprintf('Can not find %s,creating it:',SYS_config.beam_cache_file);
                fprintf('\n               ');
            end
            creat_beam = true;
        end
    else
        creat_beam = true;
    end
    % ��ʼ����������
    if creat_beam
        beam_antenna_gain=zeros(ue_len,enb_len);%��ʼ��BS�˲������;���
        beam_ue_gain=zeros(ue_len,enb_len);     %��ʼ��UE�˲������;���
        theta_map=zeros(ue_len,enb_len);        %��ʼ����BS����UE�Ĵ�ֱ��ӳ�����
        phi_map=zeros(ue_len,enb_len);          %��ʼ����BS����UE��ˮƽ��λ��ӳ�����
        theta_map_ue=zeros(ue_len,enb_len);     %��ʼ����UE����BS�Ĵ�ֱ��ӳ�����
        phi_map_ue=zeros(ue_len,enb_len);       %��ʼ����UE����BS��ˮƽ��λ��ӳ�����
        self_theta=zeros(ue_len,enb_len);       %��ʼ����BS�ڸ�UE������ʱBS�ನ��ָ��Ĵ�ֱ��
        self_phi=zeros(ue_len,enb_len);         %��ʼ����BS�ڸ�UE������ʱBS�ನ��ָ���ˮƽ��λ��
        self_phi_ue=zeros(ue_len);              %��ʼ����BS�ڸ�UE������ʱUE�ನ��ָ��Ĵ�ֱ��
        self_theta_ue=zeros(ue_len,enb_len);    %��ʼ����BS�ڸ�UE������ʱUE�ನ��ָ���ˮƽ��λ��

        for u_=1:ue_len
            if SYS_config.isDouble==false %��ϵͳ��������ϵͳ���ż���ʵ�塣����Ҫ����ֵ
                UEs(u_).attached_eNodeB.ad_interf_eNodeB_sectors.eNodeB_id=[];
            end
            
            for enb=1:enb_len
                % ���ÿ��������UE�����ȣ����˶��������л�������flag��false
                tmp= [eNodeBs_sectors(enb).attached_UEs_vector.id];
                if length(tmp)~=SYS_config.UE_per_eNodeB
                    flag=false;
                end
                if ~isempty (find([[UEs(u_).attached_eNodeB.in_interf_eNodeB_sectors.eNodeB_id],[UEs(u_).attached_eNodeB.ad_interf_eNodeB_sectors.eNodeB_id],UEs(u_).attached_eNodeB.eNodeB_id]==enb)) %���b�ڱ�ϵͳ�ڻ�վ�б������ϵͳ��վ�б��У�����Ϊ��UE�����վ
                    switch SYS_config.scene_type
                        case {'UMA', 'UMI'} % UMA��UMI BS��UEˮƽ���봹ֱ�ǵļ��㹫ʽ��ͬ
                            % UMA/UMI��ֱ�ǵļ��㣨sector(enb)��UE(u_)�Ĵ�ֱ�ǣ�
                            theta_map(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height)./pi*180;
                            % UMA/UMIˮƽ�ǵļ��㣨sector(enb)��UE(u_)��ˮƽ�ǣ�
                            phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;
                            %����ˮƽ�Ǽ����ʱ���ȥ�����߰�ķ�λ�ǣ�������Ҫ���䷶Χ����[-180��180)
                            phi_map(u_,enb)=utils.miscUtils.wrapToAll180(phi_map(u_,enb));
                            %UE(u_)��sector(enb)�ķ�λ��
                            phi_map_ue(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                            theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;
                            %��վˮƽ���漰��������UEˮƽ���漰UE��ת�Ƕȣ������Ҫ��ˮƽ�ǵķ�Χ����[-180��180)
                            phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                        case {'InH','InH2','InH3'} %InH BS�ĽǶȼ�����UMA/UMI��ͬ��UE�˵ĽǶȼ�����UMA/UMI��ͬ
                            %InH BS ��ֱ�ǵļ��㣺90-arctan(delta(y)/������ˮƽ�������վ�ߵ�ƽ����)
                            theta_map(u_,enb)=atan2(sqrt((SYS_config.site_height- SYS_config.UE_height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180;
                            %InH BS ˮƽ�ǵļ��㣺arctan(delta(y)/վ��)
                            phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),SYS_config.site_height-SYS_config.UE_height)./pi*180;                          
                            %InH UE ˮƽ���뷽��ǵļ�����UMA/UMI��ͬ��ˮƽ�Ƕ���������ת��
                            phi_map_ue(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                            theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),SYS_config.site_height-SYS_config.UE_height)./pi*180;
                            % UEˮƽ�Ǽ����ʱ���ȥ���������ת�ǣ�������Ҫ���䷶Χ����[-180��180)
                            phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                        case 'UMa_to_UMi'%UMA����UMI�����Ƕȼ��㷽����UMA/UMI��ͬ������������ϵͳ���칹�ģ���˼���Ƕ�ʹ�õĸ߶Ȳ�Ҫ��ϵͳ����
                            %BS/UEˮƽ�ǵļ�����߶Ȳ��޹أ�������Ķ�
                            phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;                            
                            phi_map(u_,enb)=utils.miscUtils.wrapToAll180(phi_map(u_,enb));
                            phi_map_ue(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                            phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                            if enb<=networkPathlossMap.num_first_sectors %��վ�ǵ�һ��ϵͳ�Ļ�վ
                                theta_map(u_,enb)=atan2(sqrt((SYS_config.site_height- UEs(u_).height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180;
                                phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;
                            else %��վ���ڵڶ���ϵͳ
                            theta_map(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height)./pi*180;
                                theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),SYS_config.site_height2-UEs(u_).height)./pi*180;
                            end
                        case{'UMi_to_InH','UMa_to_InH'} %UMA/UMI��InH�ĸ����У���һ��ϵͳʹ��InH�ĽǶȼ��㷽����ͬ��Ҫ�����칹ϵͳ�ĸ߶Ȳ������
                           %UEˮƽ�ǵļ��������������ͬ
                            phi_map_ue(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                            phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                            if enb<=networkPathlossMap.num_first_sectors %��վ�ǵ�һ��ϵͳ�Ļ�վ
                                theta_map(u_,enb)=atan2(sqrt((SYS_config.site_height- UEs(u_).height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180;
                                theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;
                                phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),SYS_config.site_height-SYS_config.UE_height)./pi*180;
                            else %��վ���ڵڶ���ϵͳ
                                theta_map(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height2)./pi*180;
                                theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),SYS_config.site_height2-UEs(u_).height)./pi*180;
                                phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;                               
                            end
                            phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                    end
                end
            end
        end
        switch SYS_config.beam_accurancy_type %�жϲ���ƫ�������
            case 'Constant' %ƫ������Ϊ����
                bs_acc_th=SYS_config.beam_accuracy_theta_bs; 
                ue_acc_th=SYS_config.beam_accuracy_theta_ue;
                bs_acc_phi=SYS_config.beam_accuracy_phi_bs;
                ue_acc_phi=SYS_config.beam_accuracy_phi_ue;
            case 'Unique' %ƫ������Ϊ���ȷֲ�
                bs_acc_th=SYS_config.beam_accuracy_theta_bs * rand(size(theta_map));
                ue_acc_th=SYS_config.beam_accuracy_theta_ue * rand(size(theta_map_ue));
                bs_acc_phi=SYS_config.beam_accuracy_phi_bs * rand(size(phi_map));
                ue_acc_phi=SYS_config.beam_accuracy_phi_ue * rand(size(phi_map_ue));
            otherwise
                error('beam accurancy type error');
        end
        %UE�˳�������ҳ���Ӧ��u_��UE��enb��Ӧ�ķ���UE�ı�ţ��ӽǶ�Map�ж�ȡ��Ҳ��ȷ����ǰ����UE������UE�Ĳ�������,��ÿ��UEָ��ǰ����С��ʱ��������С���ĽǶ�
        for u_=1:ue_len
            ser_ue_group=mod(u_,SYS_config.UE_per_eNodeB);%����UE��Ӧ�Ĳ������
            if ser_ue_group==0  %���u_���Ա�LTE_config.UE_per_eNodeB������������Ϊ�÷��񼯵����һ��С��
                ser_ue_group=SYS_config.UE_per_eNodeB;
            end
            
            for enb=1:enb_len
                if ~isempty (find([[UEs(u_).attached_eNodeB.in_interf_eNodeB_sectors.eNodeB_id],[UEs(u_).attached_eNodeB.ad_interf_eNodeB_sectors.eNodeB_id],UEs(u_).attached_eNodeB.eNodeB_id]==enb)) %���b�ڱ�ϵͳ�ڻ�վ�б������ϵͳ��վ�б��У�����Ϊ��UE�����վ
                    %��idΪu_��UE������ʱ,UE(u_)�Ĳ���ָ�������С��
                    self_theta_ue(u_,enb)=theta_map_ue(u_,UEs(u_).attached_eNodeB.eNodeB_id);
                    self_phi_ue(u_,enb)=phi_map_ue(u_,UEs(u_).attached_eNodeB.eNodeB_id);
                    if flag %���ÿ��С�����û���һ��������÷������ʽ����slot
                        %��idΪu_��UE������ʱ,sector(enb)�Ĳ���ָ�������ue
                        self_theta(u_,enb)=theta_map(ser_ue_group+(enb-1)*SYS_config.UE_per_eNodeB,enb);
                        %�ҳ���Ӧ��u_��UE��enb��Ӧ�ķ���UE�ı�ţ��ӽǶ�Map�ж�ȡ��Ҳ��ȷ����ǰ����UE������sector�Ĳ�������,��ÿ��SECTORָ���伯�����뵱ǰUE�����ͬ�ı���UE
                        self_phi(u_,enb)=phi_map(ser_ue_group+(enb-1)*SYS_config.UE_per_eNodeB,enb);
                    else %���С�����û�����һ���������С������������û��ķ�ʽ
                        tmp= [eNodeBs_sectors(enb).attached_UEs_vector.id]; %enbС���ķ���UE ��id��
                        serve_ue_id=tmp(randi(length(tmp))); %�ӷ��������ѡȡһ��UE                     
                        self_theta(u_,enb)=theta_map(serve_ue_id,enb);
                        %�ҳ���Ӧ��u_��UE��enb��Ӧ�ķ���UE�ı�ţ��ӽǶ�Map�ж�ȡ��Ҳ��ȷ����ǰ����UE������sector�Ĳ�������,��ÿ��SECTORָ���伯�����뵱ǰUE�����ͬ�ı���UE
                        self_phi(u_,enb)=phi_map(serve_ue_id,enb);
                    end
                    
                end
            end
        end
        %% 4G����36942��ģ�ͣ������ǲ������͡�NR����38900��38901��ģ�ͣ����ò�������
        % NR����LTE
        if strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && ~strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')%�����һ��ϵͳ��LTE
            for u_ = 1:ue_len %��վ�ನ�����͵ļ��㣬LTE��վ������UE�Ĳ�����������Ϊ0��NR��������
                if ~isempty(find(u_==u_markings,1)) %���ý�����
                    if DEBUG_LEVEL>=1
                        fprintf('=');
                        if u_ == ue_len;
                            fprintf('\n');
                        end
                    end
                end
                for enb = 1:networkPathlossMap.num_first_sectors %��һ��ϵͳ��վ������UE�Ļ�վ�ನ������Ϊ0
                    beam_antenna_gain(u_,enb)=0;
                end
                for enb = networkPathlossMap.num_first_sectors+1:enb_len
                    if SYS_config.antenna_mode==1||SYS_config.antenna_mode==2 %���BS���ò�������
                        beam_antenna_gain(u_,enb)=UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_map(u_,enb),phi_map(u_,enb),self_theta(u_,enb)-90+bs_acc_th,self_phi(u_,enb)+bs_acc_phi)+3;%����3dB��������
                    end
                end
            end
            for enb = 1:enb_len %UE�ನ�����͵ļ��㣬lteUE������UE�Ĳ�����������Ϊ0��NR��UE��������
                for u_ = 1:networkPathlossMap.num_first_UEs %��һ��ϵͳUE������С����UE�ನ������Ϊ0
                    beam_ue_gain(u_,enb)=0;
                end
                for u_ = networkPathlossMap.num_first_UEs+1:ue_len
                    if SYS_config.antenna_mode==2 %���UE���ò�������
                        beam_ue_gain(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_map_ue(u_,enb),phi_map_ue(u_,enb),self_theta_ue(u_,enb)-90+ue_acc_th,self_phi_ue(u_,enb)+ue_acc_phi);
                    end
                end
            end
            %LTE����NR
        elseif  ~strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')%����ڶ���ϵͳ��LTE
            for u_ = 1:ue_len %��վ�ನ�����͵ļ��㣬LTE��վ������UE�Ĳ�����������Ϊ0��NR��������
                if ~isempty(find(u_==u_markings,1)) %���ý�����
                    if DEBUG_LEVEL>=1
                        fprintf('=');
                        if u_ == ue_len;
                            fprintf('\n');
                        end
                    end
                end
                for enb = 1:networkPathlossMap.num_first_sectors %������һ��С���Ļ�վ
                    if SYS_config.antenna_mode==1||SYS_config.antenna_mode==2 %���BS���ò�������
                        beam_antenna_gain(u_,enb)=UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_map(u_,enb),phi_map(u_,enb),self_theta(u_,enb)-90+bs_acc_th,self_phi(u_,enb)+bs_acc_phi)+3;%����3dB��������
                    end
                end
                for enb = networkPathlossMap.num_first_sectors+1:enb_len %�ڶ���ϵͳ��վ������UE�Ļ�վ�ನ������Ϊ0
                    beam_antenna_gain(u_,enb)=0;
                end
            end
            for enb = 1:enb_len %UE�ನ�����͵ļ��㣬lteUE������UE�Ĳ�����������Ϊ0��NR��UE��������
                for u_ = 1:networkPathlossMap.num_first_UEs %������һ��С��UE
                    if SYS_config.antenna_mode==2
                        beam_ue_gain(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_map_ue(u_,enb),phi_map_ue(u_,enb),self_theta_ue(u_,enb)-90+ue_acc_th,self_phi_ue(u_,enb)+ue_acc_phi);
                    end
                    
                end
                for u_ = networkPathlossMap.num_first_UEs+1:ue_len %��һ��ϵͳUE������С����UE�ನ������Ϊ0
                    beam_ue_gain(u_,enb)=0;
                end
            end
            %LTE����LTE
        elseif strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')
            if DEBUG_LEVEL>=1
                fprintf('non NR, no beam forming\n');
            end
            %NR ����NR
        elseif ~strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && ~strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')
            for u_=1:ue_len
                for enb=1:enb_len
                    if SYS_config.antenna_mode==1||SYS_config.antenna_mode==2 %���BS���ò�������
                        beam_antenna_gain(u_,enb)=UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_map(u_,enb),phi_map(u_,enb),self_theta(u_,enb)-90+bs_acc_th,self_phi(u_,enb)+bs_acc_phi)+3;%����3dB��������
                    end
                    if SYS_config.antenna_mode==2 %���UE���ò�������
                        beam_ue_gain(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_map_ue(u_,enb),phi_map_ue(u_,enb),self_theta_ue(u_,enb)-90+ue_acc_th,self_phi_ue(u_,enb)+ue_acc_phi);
                    end
                end
                if ~isempty(find(u_==u_markings,1)) %���ý�����
                    if DEBUG_LEVEL>=1
                        fprintf('=');
                        if u_ == ue_len;
                            fprintf('\n');
                        end
                    end
                end
            end
        end
        %���沨�����;���
        if SYS_config.isSave
            filename1 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain.mat'];
            filename2 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain.mat'];
            save(filename1,'beam_antenna_gain');
            save(filename2,'beam_ue_gain');
        end
    end
end
end

