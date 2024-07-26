function [ue_or_bs_gain_asy,ue_gain_asy,bs_or_ue_gain_asy,bs_gain_asy] = NR_beam_forming_asy(UEs,eNodeBs_sectors,networkPathlossMap,SYS_config)
% ���㲻ͬ�������µĲ������;���һ����Ҫ˫ϵͳ��
% ������������ͼ��
% ���������
% UEs���û��б�
% eNodeBs_sectors��С��ʵ��
% networkPathlossMap��·��ӳ�����
% SYS_config��LTE����ϵͳ
%
% ���������
% ue_or_bs_gain_asy�����и��������У�����˲������;��󡪡���ϵͳϵͳ��BS�˲����������棬��ϵͳ�ж�ӦΪ����UE������������
% ue_gain_asy�����и��������У����ն˲������;��󡪡�UE�˲�����������
% bs_or_ue_gain_asy�����и��������У�����˲������;��󡪡���ϵͳ��UE�����������棬��ϵͳ�ж�ӦΪ����BS������������
% bs_gain_asy�� ���и��������У����ն˲������;��󡪡�BS�˲�����������
flag=true; %���ڱ�ʶ��ϵͳÿBS��UE���Ƿ���ͬ
DEBUG_LEVEL = SYS_config.debug_level;
num_markers = 10; %����������
u_markings  = round(linspace(1,length(UEs),num_markers));
ue_len=length(UEs);
enb_len=length(eNodeBs_sectors);%����ϵͳ�Ļ������122����վ������UMA wrap around��,������sector����ֵΪ366
ue_or_bs_gain_asy=zeros(ue_len,enb_len );
ue_gain_asy=zeros(ue_len,enb_len );
bs_or_ue_gain_asy=zeros(ue_len,enb_len );
bs_gain_asy=zeros(ue_len,enb_len );

if DEBUG_LEVEL>=1
    if SYS_config.antenna_mode==2
        fprintf('Creating BS and UE beam forming for asynchronization:');
    elseif SYS_config.antenna_mode==1
        fprintf('Creating BS beam forming for asynchronization:');
    end
end


if SYS_config.antenna_mode==0
    if DEBUG_LEVEL>=1
        fprintf('No BF in BS nor UE\n')
    end
else
    %���벨����������
    if SYS_config.use_cache
        file_name1 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain_up_down.mat'];
        file_name2 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain_up_down.mat'];
        file_name3 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain_down_up.mat'];
        file_name4 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain_down_up.mat'];
        try
            tmp = load(file_name1);
            ue_or_bs_gain_asy = tmp.ue_or_bs_gain_asy;
            tmp = load(file_name2);
            ue_gain_asy = tmp.ue_gain_asy;
            tmp = load(file_name3);
            bs_gain_asy = tmp.bs_gain_asy;
            tmp = load(file_name4);
            bs_or_ue_gain_asy = tmp.bs_or_ue_gain_asy;
            creat_beam = false;
            if DEBUG_LEVEL>=1
                fprintf('(use beam forming cache)\n');
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
    %% ��ʼ������
    if creat_beam
        ue_or_bs_gain_asy=zeros(ue_len,enb_len);     %��ʼ�����и�������ʱ����˲������ξ��󣨱�ϵͳBS�����źţ���ϵͳUE�����źţ�
        ue_gain_asy=zeros(ue_len,enb_len);           %��ʼ�����и�������ʱ���ն˲������;��󣨱�ϵͳ����ϵͳUE�˽����źţ�
        bs_or_ue_gain_asy=zeros(ue_len,enb_len);     %��ʼ�����и�������ʱ����˲������;��󣨱�ϵͳUE�����źţ���ϵͳBS�����źţ�
        bs_gain_asy=zeros(ue_len,enb_len);           %��ʼ�����и�������ʱ���ն˲������ξ��󣨱�ϵͳ����ϵͳBS�˽����źţ�
        theta_map_ue=zeros(ue_len,enb_len);          %��ʼ����UE��������ϵͳ����վ�Ĵ�ֱ�Ǿ���ϵͳ�ڸ�����ʽ���䣩
        phi_map_ue=zeros(ue_len,enb_len);            %��ʼ����UE��������ϵͳ����վ��ˮƽ�Ǿ���
        theta_map=zeros(ue_len,enb_len);             %��ʼ����BS��������ϵͳ��UE�Ĵ�ֱ�Ǿ���
        phi_map=zeros(ue_len,enb_len);               %��ʼ����BS��������ϵͳ��UE��ˮƽ�Ǿ���
        self_theta=zeros(ue_len,enb_len);
        self_theta_ue=zeros(ue_len,enb_len);
        self_phi=zeros(ue_len,enb_len);
        self_phi_ue=zeros(ue_len);
        beam_orientation_phi_ue2bs=zeros(ue_len);    %��ʼ��UE�������BS�Ĳ���ˮƽ�����
        beam_orientation_phi_bs2ue=zeros(ue_len);    %��ʼ������BS��UE�Ĳ���ˮƽ�����
        s=zeros(ue_len);                             %UE�����վ��id
        beam_orientation_theta_ue2bs=zeros(ue_len);  %��ʼ��UE�������BS�Ĳ�������ֱ��
        beam_orientation_theta_bs2ue=zeros(ue_len);  %��ʼ������BS��UE�Ĳ�����ֱ��
        ser_group_num=zeros(ue_len);                 %�������
        h1=zeros(ue_len);                            %��ʼ��vicϵͳ��վ�ĸ߶�
        h2=zeros(ue_len);                            %��ʼ��aggϵͳ��վ�ĸ߶�
        phi_up_down_agg=zeros(ue_len,enb_len);       %��ʼ��aggϵͳUE��vicϵͳUE��ˮƽ��
        phi_up_down_vic=zeros(ue_len,enb_len);       %��ʼ��vicϵͳUE��aggϵͳUE��ˮƽ��
        theta_up_down_agg=zeros(ue_len,enb_len);     %��ʼ��aggϵͳUE��vicϵͳUE�Ĵ�ֱ��
        theta_up_down_vic=zeros(ue_len,enb_len);     %��ʼ��vicϵͳUE��aggϵͳUE�Ĵ�ֱ��
        phi_down_up_agg=zeros(ue_len,enb_len);       %��ʼ��aggϵͳBS��vicϵͳBS��ˮƽ��
        phi_down_up_vic=zeros(ue_len,enb_len);       %��ʼ��vicϵͳBS��aggϵͳBS��ˮƽ��
        theta_down_up_agg=zeros(ue_len,enb_len);     %��ʼ��aggϵͳBS��vicϵͳBS�Ĵ�ֱ��
        theta_down_up_vic=zeros(ue_len,enb_len);     %��ʼ��vicϵͳBS��aggϵͳBS�Ĵ�ֱ��
        %��ÿ��UE����UE�����ߣ�������ÿ��UE�Ĳ�������
        for enb=1:enb_len
            % ���ÿ��������UE�����ȣ����˶��������л�������flag��false
            tmp= [eNodeBs_sectors(enb).attached_UEs_vector.id];
            if length(tmp)~=SYS_config.UE_per_eNodeB
                flag=false;
            end
        end
        
        for u_=1:ue_len
            s(u_)=UEs(u_).attached_eNodeB.eNodeB_id; %����ȷ����UE�����sector
            beam_orientation_phi_ue2bs(u_)=atan2(eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;%����UEs(u_)��������վ��ˮƽ��
            %�������ΪInH����ΪUMa_to_InH��UMi_to_InH�����µ�һ��ϵͳʱ������InH�ĽǶȼ��㷽��
            if strcmp ('InH',SYS_config.scene_type) || strcmp ('InH2',SYS_config.scene_type) ||strcmp ('InH3',SYS_config.scene_type) ||( u_<= networkPathlossMap.num_first_UEs && (strcmp('UMa_to_InH',SYS_config.scene_type) || strcmp('UMi_to_InH',SYS_config.scene_type)))
                beam_orientation_phi_bs2ue(u_)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2),SYS_config.site_height-SYS_config.UE_height)./pi*180;
            else
                beam_orientation_phi_bs2ue(u_)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1))./pi*180-UEs(u_).orientation;%����������վ��UE(u_)��ˮƽ��
            end
            % ������ˮƽ��ת��Ϊ(-180,180]
            beam_orientation_phi_ue2bs(u_)=utils.miscUtils.wrapToAll180(beam_orientation_phi_ue2bs(u_));
            beam_orientation_phi_bs2ue(u_)=utils.miscUtils.wrapToAll180(beam_orientation_phi_bs2ue(u_));
            switch SYS_config.scene_type
                case {'InH','InH2','InH3'}
                    beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;%����UE(u_)��������վ�Ĵ�ֱ��
                    beam_orientation_theta_bs2ue(u_)=atan2(sqrt((SYS_config.site_height- UEs(u_).height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1))./pi*180;%��������վ��UE(u_)�Ĵ�ֱ��
                case {'UMa_to_InH','UMi_to_InH'}
                    if u_<= networkPathlossMap.num_first_UEs %����ǵ�һ��ϵͳ������InH�����������Ƕ�
                        beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;%����UE(u_)��������վ�Ĵ�ֱ��
                        beam_orientation_theta_bs2ue(u_)=atan2(sqrt((SYS_config.site_height- UEs(u_).height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1))./pi*180;%��������վ��UE(u_)�Ĵ�ֱ��
                    else %����ǵڶ���ϵͳ
                        beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height2-UEs(u_).height)./pi*180;%����UE(u_)��������վ�Ĵ�ֱ��
                        beam_orientation_theta_bs2ue(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height2)./pi*180;%����UE(u_)��������վ�Ĵ�ֱ��
                    end
                case 'UMa_to_UMi'
                    if u_<= networkPathlossMap.num_first_UEs %UE���Ե�һ��ϵͳ
                        beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;%����UE(u_)��������վ�Ĵ�ֱ��
                        beam_orientation_theta_bs2ue(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height)./pi*180;%����UE(u_)��������վ�Ĵ�ֱ��
                    else %UE���Եڶ���ϵͳ
                        beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height2-UEs(u_).height)./pi*180;%����UE(u_)��������վ�Ĵ�ֱ��
                        beam_orientation_theta_bs2ue(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height2)./pi*180;%����UE(u_)��������վ�Ĵ�ֱ��
                    end
                otherwise
                    beam_orientation_theta_ue2bs(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),SYS_config.site_height-UEs(u_).height)./pi*180;%����UE(u_)��������վ�Ĵ�ֱ��
                    beam_orientation_theta_bs2ue(u_)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(s(u_)).parent_eNodeB.pos(2)).^2),UEs(u_).height-SYS_config.site_height)./pi*180;%����UE(u_)��������վ�Ĵ�ֱ��

            end
            
            if flag % ���ÿ��sector���û�����ͬ�������UEs(u_)�ķ������
                ser_group_num(u_)=mod(u_,SYS_config.UE_per_eNodeB);
                if ser_group_num(u_)==0 %���ȡ��Ϊ0�����������ΪLTE_config.UE_per_eNodeB
                    ser_group_num(u_)=SYS_config.UE_per_eNodeB;
                end
            else %���ÿ��sector�û�����ͬ���������ÿ���û���UE�б���ѡ��һ��
                tmp= [UEs(u_).attached_eNodeB.attached_UEs_vector.id];
                rand_index=randi(length(tmp));
                serve_id(u_)=tmp(rand_index);
            end
            
            
        end
        %����UE�������һ��UE�Ĳ�������
        for u_=1:ue_len
            for enb=1:enb_len
                if (u_<= networkPathlossMap.num_first_UEs &&enb>networkPathlossMap.num_first_sectors) || (u_>networkPathlossMap.num_first_UEs && enb<=networkPathlossMap.num_first_sectors)
                    if enb>networkPathlossMap.num_first_sectors %��վ���Եڶ���ϵͳ
                       h1(u_)=SYS_config.site_height; %h1��ʾvicϵͳ��վ�ĸ߶�
                       switch SYS_config.scene_type
                           case {'UMa_to_UMi','UMa_to_InH', 'UMi_to_InH'}
                               h2(u_)=SYS_config.site_height2; %h2��ʾaggϵͳ��վ�߶�
                           otherwise
                               h2(u_)=SYS_config.site_height;
                       end
                   else %���Ż�վ���Ե�һ��ϵͳ���ڶ���ϵͳ��vicϵͳ
                       h2(u_)=SYS_config.site_height;
                       switch SYS_config.scene_type
                           case {'UMa_to_UMi','UMa_to_InH', 'UMi_to_InH'}
                               h1(u_)=SYS_config.site_height2;
                           otherwise
                               h1(u_)=SYS_config.site_height;
                       end
                   end
                   
                if flag                                       
                    %����aggϵͳUE��vicϵͳUE��ˮƽ��(���и�������,�Ӹ���UE�ķ����������)
                    phi_up_down_agg(u_,enb)=atan2(UEs(u_).pos(2)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(2),UEs(u_).pos(1)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(1))./pi*180-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).orientation;%����UE��UE(u_)���丳������
                    %����vicϵͳUE��aggϵͳUE��ˮƽ��(���и�������,�ӱ�����UE�Ľ����������)
                    phi_up_down_vic(u_,enb)=atan2(UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(2)-UEs(u_).pos(2),UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                    %����aggϵͳUE��vicϵͳUE�Ĵ�ֱ��
                    theta_up_down_agg(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(1)).^2+(UEs(u_).pos(2)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(2)).^2),UEs(u_).height-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).height)./pi*180;
                    %����vicϵͳUE��aggϵͳUE�Ĵ�ֱ��
                    theta_up_down_vic(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(1)).^2+(UEs(u_).pos(2)-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).pos(2)).^2),UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).height-UEs(u_).height)./pi*180;
                else 
                    phi_up_down_agg(u_,enb)=atan2(UEs(u_).pos(2)-UEs(serve_id(u_)).pos(2),UEs(u_).pos(1)-UEs(serve_id(u_)).pos(1))./pi*180-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).orientation;%����UE��UE(u_)���丳������
                    phi_up_down_vic(u_,enb)=atan2(UEs(serve_id(u_)).pos(2)-UEs(u_).pos(2),UEs(serve_id(u_)).pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                    theta_up_down_agg(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-UEs(serve_id(u_)).pos(1)).^2+(UEs(u_).pos(2)-UEs(serve_id(u_)).pos(2)).^2),UEs(u_).height-UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).height)./pi*180;
                    theta_up_down_vic(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-UEs(serve_id(u_)).pos(1)).^2+(UEs(u_).pos(2)-UEs(serve_id(u_)).pos(2)).^2),UEs(serve_id).height-UEs(u_).height)./pi*180;
                end

                    %����aggϵͳsector��vicϵͳsector��ˮƽ�ǣ����и�������,�Ӹ���sector�ķ������������
                    phi_down_up_agg(u_,enb)=atan2(UEs(u_).attached_eNodeB.parent_eNodeB.pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).attached_eNodeB.parent_eNodeB.pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))-eNodeBs_sectors(enb).azimuth;
                    %����vicϵͳsector��aggϵͳsector��ˮƽ�ǣ����и�������,�ӱ�����sector�Ľ������������
                    phi_down_up_vic(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).attached_eNodeB.parent_eNodeB.pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).attached_eNodeB.parent_eNodeB.pos(1))-UEs(u_).attached_eNodeB.azimuth;
                    %�Ƕ�ת��Ϊ(-180,180]
                    phi_up_down_agg(u_,enb)=utils.miscUtils.wrapToAll180(phi_up_down_agg(u_,enb));
                    phi_up_down_vic(u_,enb)=utils.miscUtils.wrapToAll180(phi_up_down_vic(u_,enb));
                    phi_down_up_agg(u_,enb)=utils.miscUtils.wrapToAll180(phi_down_up_agg(u_,enb));
                    phi_down_up_vic(u_,enb)=utils.miscUtils.wrapToAll180(phi_down_up_vic(u_,enb));
    
                    %����aggϵͳsector��vicϵͳsector�Ĵ�ֱ��
                    theta_down_up_agg(u_,enb)=atan2(sqrt((UEs(u_).attached_eNodeB.parent_eNodeB.pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).attached_eNodeB.parent_eNodeB.pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),h1(u_)-h2(u_))./pi*180;
                    %����vicϵͳsector��vicϵͳsector�Ĵ�ֱ��
                    theta_down_up_vic(u_,enb)=atan2(sqrt((eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).attached_eNodeB.parent_eNodeB.pos(1)).^2+(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).attached_eNodeB.parent_eNodeB.pos(2)).^2),h2(u_)-h1(u_))./pi*180;

                else
                    %UE�˿���ˮƽ��ת�ǣ����Ҵ�ֱ����ˮƽ�����վ�˵Ļ���
                    phi_map_ue(u_,enb)=atan2(eNodeBs_sectors(enb).parent_eNodeB.pos(2)-UEs(u_).pos(2),eNodeBs_sectors(enb).parent_eNodeB.pos(1)-UEs(u_).pos(1))./pi*180-UEs(u_).orientation;
                    theta_map_ue(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),h1(u_)-UEs(u_).height)./pi*180;
                    
                     %h1��¼����victim��ϵͳվ�ߣ�������������һ��ϵͳ��UEʱ����һ��ϵͳΪvictimϵͳ���������ڶ���ϵͳ��UEʱ���ڶ���ϵͳΪ��Ӧ��"victimϵͳ"
                     %�����൱����h1Ϊu_����ϵͳ��վ�ߡ�
                    switch SYS_config.scene_type
                        case {'UMA', 'UMI','RMa','UMa_to_UMi'}
                            % UMA/UMI��ֱ�ǵļ���
                            theta_map(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).height-h1(u_))./pi*180;
                            % UMA/UMIˮƽ�ǵļ��㣺arctan(delta(y)/delta(x))-sector azimuth
                            phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;
                        case {'InH','InH2','InH3'}
                            %InH BS ��ֱ�ǵļ��㣺90-arctan(delta(y)/������ˮƽ�������վ�ߵ�ƽ����)
                            theta_map(u_,enb)=atan2(sqrt((SYS_config.site_height- SYS_config.UE_height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180;
                            %InH BS ˮƽ�ǵļ��㣺arctan(delta(y)/վ��)
                            phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),h1(u_)-SYS_config.UE_height)./pi*180;                            
                        case {'UMa_to_InH','UMi_to_InH'}
                            if enb <= networkPathlossMap.num_first_sectors %��һ��ϵͳ����InH�ķ�ʽ
                                theta_map(u_,enb)=atan2(sqrt((h1(u_)- SYS_config.UE_height).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180;
                                phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),h1(u_)-UEs(u_).height)./pi*180;
                            else %�ڶ���ϵͳ����UMA��ʽ
                                theta_map(u_,enb)=atan2(sqrt((UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1)).^2+(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2)).^2),UEs(u_).height-h1(u_))./pi*180;
                                phi_map(u_,enb)=atan2(UEs(u_).pos(2)-eNodeBs_sectors(enb).parent_eNodeB.pos(2),UEs(u_).pos(1)-eNodeBs_sectors(enb).parent_eNodeB.pos(1))./pi*180-eNodeBs_sectors(enb).azimuth;
                            end
                                
                    end
                    % UEˮƽ�Ǽ����ʱ���ȥ���������ת�ǣ�������Ҫ���䷶Χ����[-180��180)
                    phi_map_ue(u_,enb)=utils.miscUtils.wrapToAll180(phi_map_ue(u_,enb));
                    %����ˮƽ�Ǽ����ʱ���ȥ�����߰�ķ�λ�ǣ�������Ҫ���䷶Χ����[-180��180)
                    phi_map(u_,enb)=utils.miscUtils.wrapToAll180(phi_map(u_,enb));
                end
            end
        end
        
        switch SYS_config.beam_accurancy_type %�жϲ���ƫ�������
            case 'Constant' %ƫ������Ϊ����
                bs_acc_th=SYS_config.beam_accuracy_theta_bs;
                ue_acc_th=SYS_config.beam_accuracy_theta_ue;
                bs_acc_phi=SYS_config.beam_accuracy_phi_bs;
                ue_acc_phi=SYS_config.beam_accuracy_phi_ue;
            case 'Unique'  %ƫ������Ϊ���ȷֲ�
                bs_acc_th=SYS_config.beam_accuracy_theta_bs * rand(size(theta_map));
                ue_acc_th=SYS_config.beam_accuracy_theta_ue * rand(size(theta_map_ue));
                bs_acc_phi=SYS_config.beam_accuracy_phi_bs * rand(size(phi_map));
                ue_acc_phi=SYS_config.beam_accuracy_phi_ue * rand(size(phi_map_ue));
            otherwise
                error('beam accurancy type error')
        end
        
        for u_=1:ue_len
            for enb=1:enb_len
                if (u_<= networkPathlossMap.num_first_UEs &&enb>networkPathlossMap.num_first_sectors) || (u_>networkPathlossMap.num_first_UEs && enb<=networkPathlossMap.num_first_sectors)
                    if SYS_config.antenna_mode==2
                        %�����и������У�UE(u_)�Ľ��ո������棨�����ϵͳ��
                        ue_gain_asy(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_up_down_vic(u_,enb),phi_up_down_vic(u_,enb),beam_orientation_theta_ue2bs(u_)-90+ue_acc_th,beam_orientation_phi_ue2bs(u_)+ue_acc_phi);
                        %�����и������У�UE(u_)������sector���ո������棨�����ϵͳ��
                        bs_gain_asy(u_,enb) = UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_down_up_vic(u_,enb),phi_down_up_vic(u_,enb),beam_orientation_theta_bs2ue(u_)-90+bs_acc_th,beam_orientation_phi_bs2ue(u_)+bs_acc_phi);
                        if flag
                            %�����и������У�UE(u_)�ܵ���UE����ϵͳ�����ŷ��丳������
                            ue_or_bs_gain_asy(u_,enb)=UEs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_)).antenna.beam_gain(1,theta_up_down_agg((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_),s(u_)),phi_up_down_agg((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_),s(u_)),beam_orientation_theta_ue2bs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_))-90+ue_acc_th,beam_orientation_phi_ue2bs((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_))+ue_acc_phi);
                            %�����и������У�UE(u_)������sector�ܵ�sector����ϵͳ�����ŷ��丳������
                            bs_or_ue_gain_asy(u_,enb) = eNodeBs_sectors(enb).antenna.beam_gain(1,theta_down_up_agg((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_),s(u_)),phi_down_up_agg((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_),s(u_)),beam_orientation_theta_bs2ue((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_))-90+bs_acc_th,beam_orientation_phi_bs2ue((enb-1)*SYS_config.UE_per_eNodeB+ser_group_num(u_))+bs_acc_phi);
                        else
                            ue_or_bs_gain_asy(u_,enb)=UEs(serve(u_)).antenna.beam_gain(1,theta_up_down_agg(serve_id(u_),s(u_)),phi_up_down_agg(serve_id(u_),s(u_)),beam_orientation_theta_ue2bs(serve_id(u_))-90+ue_acc_th,beam_orientation_phi_ue2bs(serve_id(u_))+ue_acc_phi);
                            bs_or_ue_gain_asy(u_,enb) = eNodeBs_sectors(enb).antenna.beam_gain(1,theta_down_up_agg(serve_id(u_),s(u_)),phi_down_up_agg(serve_id(u_),s(u_)),beam_orientation_theta_bs2ue(serve_id(u_))-90+bs_acc_th,beam_orientation_phi_bs2ue(serve_id(u_))+bs_acc_phi);
                        end
                    end
                 else
                    % ���㵱ǰUE��u_��������ʱ��enb�������UE���ɵĽ�
                    if ~flag %�����ϵͳÿ����վ������û�����һ���������ѡ�����UEȡ����
                        self_theta(u_,enb)=theta_map(serve_id(u_),enb);
                        self_phi(u_,enb)=phi_map(serve_id(u_),enb);
                    else %�����ϵͳÿ����վ������û���һ��������ѡ�����UEȡ����
                        self_theta(u_,enb)=theta_map(ser_group_num(u_)+(enb-1)*SYS_config.UE_per_eNodeB,enb);
                        self_phi(u_,enb)=phi_map(ser_group_num(u_)+(enb-1)*SYS_config.UE_per_eNodeB,enb);
                    end
                    self_theta_ue(u_,enb)=theta_map_ue(u_,s(u_));
                    self_phi_ue(u_,enb)=phi_map_ue(u_,s(u_));
                    if SYS_config.antenna_mode==2
                        ue_gain_asy(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_map_ue(u_,enb),phi_map_ue(u_,enb),self_theta_ue(u_,enb)-90+ue_acc_th,self_phi_ue(u_,enb)+ue_acc_phi);
                        bs_or_ue_gain_asy(u_,enb)=UEs(u_).antenna.beam_gain(1,theta_map_ue(u_,enb),phi_map_ue(u_,enb),self_theta_ue(u_,enb)-90+ue_acc_th,self_phi_ue(u_,enb)+ue_acc_phi);
                    end
                    if SYS_config.antenna_mode==1||SYS_config.antenna_mode==2
                        ue_or_bs_gain_asy(u_,enb)= UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_map(u_,enb),phi_map(u_,enb),self_theta(u_,enb)-90+bs_acc_th,self_phi(u_,enb)+bs_acc_phi)+3;%����3dB��������
                        bs_gain_asy(u_,enb)= UEs(u_).attached_eNodeB.antenna.beam_gain(1,theta_map(u_,enb),phi_map(u_,enb),self_theta(u_,enb)-90+bs_acc_th,self_phi(u_,enb)+bs_acc_phi)+3;%����3dB��������
                    end
                end
            end
            if ~isempty(find(u_==u_markings,1)) %���ý�����
                if DEBUG_LEVEL>=1
                    fprintf('=');
                    if u_ == ue_len
                        fprintf('\n');
                    end
                end
            end
        end
    end
    %NR����LTE
    if strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && ~strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')%�����һ��ϵͳ��LTE
        for u_=1:ue_len
            if ~isempty(find(u_==u_markings,1)) %���ý�����
                if DEBUG_LEVEL>=1
                    fprintf('=');
                    if u_ == ue_len
                        fprintf('\n');
                    end
                end
            end
            for enb = 1:networkPathlossMap.num_first_sectors   %��һ��ϵͳ��վ�Ե�һ��ϵͳ��UE�Ļ�վ�ำ������Ϊ0����һ��ϵͳUE�Եڶ���ϵͳ��UE�Ľ��ո�������Ϊ0
                ue_or_bs_gain_asy(u_,enb)=0;
                bs_or_ue_gain_asy(u_,enb)=0;
            end
        end
        for enb=1:enb_len
            for u_=1:networkPathlossMap.num_first_UEs %��һ��ϵͳUE�Ե�һ��ϵͳ�Ļ�վ��UE�ำ������Ϊ0����һ��ϵͳUE�Եڶ���ϵͳ�ķ��丳������Ϊ0
                ue_gain_asy(u_,enb)=0;
                bs_gain_asy(u_,enb)=0;
            end
        end
        %LTE����NR
    elseif ~strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')%����ڶ���ϵͳ��LTE
        for u_=1:ue_len
            if ~isempty(find(u_==u_markings,1)) %���ý�����
                if DEBUG_LEVEL>=1
                    fprintf('=');
                    if u_ == ue_len
                        fprintf('\n');
                    end
                end
            end
            for enb = networkPathlossMap.num_first_sectors+1:enb_len   %�ڶ���ϵͳ��վ�Եڶ���ϵͳ��UE�Ļ�վ�ำ������Ϊ0���ڶ���ϵͳUE�Ե�һ��ϵͳ��UE�Ľ��ո�������Ϊ0
                ue_or_bs_gain_asy(u_,enb)=0;
                bs_or_ue_gain_asy(u_,enb)=0;
            end
        end
        for enb=1:enb_len
            for u_=networkPathlossMap.num_first_UEs+1:ue_len    %�ڶ���ϵͳUE�Եڶ���ϵͳ�Ļ�վ��UE�ำ������Ϊ0���ڶ���ϵͳUE�Ե�һ��ϵͳ�ķ��丳������Ϊ0
                ue_gain_asy(u_,enb)=0;
                bs_gain_asy(u_,enb)=0;
            end
        end
    elseif strcmp(SYS_config.macroscopic_pathloss_model,'TS36942') && strcmp(SYS_config.macroscopic_pathloss_model2,'TS36942')
        if DEBUG_LEVEL>=1
            fprintf('non NR, no beam forming\n');
        end
        ue_or_bs_gain_asy=zeros(ue_len,enb_len);
        ue_gain_asy=zeros(ue_len,enb_len);
    end
    
    %���沨�����;���
    if SYS_config.isSave
        file_name1 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain_up_down.mat'];
        file_name2 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain_up_down.mat'];
        file_name3 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamAntennaGain_down_up.mat'];
        file_name4 = [SYS_config.beam_cache_folder,SYS_config.beam_cache_file,'_BeamUeGain_down_up.mat'];
        save(file_name1,'ue_or_bs_gain_asy');
        save(file_name2,'ue_gain_asy');
        save(file_name3,'bs_gain_asy');
        save(file_name4,'bs_or_ue_gain_asy');
    end
end

end
