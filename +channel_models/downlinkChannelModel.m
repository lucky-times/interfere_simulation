classdef downlinkChannelModel < handle
% ָ��UE�������ŵ�ģ�͡���ÿ��UE����һ���ض����ŵ�ģ�ͣ�����һЩ��·��Ϣ��
% ֮���Խ�downlinkChannelModel������Ϊ·����Ӱ˥�䡢��͸��ĵ������Ϣ������������·��˵����ȵ�
% ����Ķ�����ΪUE��һ�����Ա�����UE�У���һ��UE����һ����·��Ϣ

   properties

       macroscopic_pathloss_model_is_set = false;
       macroscopic_pathloss_model % ·��ģ��
       shadow_fading_model_is_set = false; 
       shadow_fading_model % ��Ӱ˥��ģ�� 
       attached_UE % ���ŵ�������UE
   end

   methods
       %% ���캯��
       function obj = downlinkChannelModel(aUE)
           obj.attached_UE = aUE;
       end
       
       %% ����һϵ�к����������ǵ���macroscopicPathlossMap��shadowFadingMapClaussen�ĺ�����ȡ�������
       % ���巽���ǣ�ͨ��attached_UE����UE�е�id��pos��������С����������Ϊ��������macroscopicPathlossMap��shadowFadingMapClaussen
       % �еĲ�������Щ����ž�����ĵĺ�����Ҫ����һϵ�е�UE����interferingEnodeBids��һ��1*Nά�ľ��󣩣���������񾶵ķֿ�д
       
       %% �õ�UE�ĸ߶�
         function height = UE_height(obj)
           pos = obj.attached_UE.pos;
           height = obj.macroscopic_pathloss_model.get_UE_height(pos);
         end
         %% �õ�UE����ת��
         function orientation = UE_orientation(obj)
           pos = obj.attached_UE.pos;
           orientation = obj.macroscopic_pathloss_model.get_UE_orientation(pos);
         end
        %% �õ�UE�������С���Ĵ�͸��� 
       function penetration_loss = macroscopic_penetration_loss(obj)
           eNodeB_id = obj.attached_UE.attached_eNodeB.eNodeB_id;
           pos = obj.attached_UE.pos;
           penetration_loss = obj.macroscopic_pathloss_model.get_penetration_loss_eNodeB(pos,eNodeB_id);
       end
       %% �õ�UE��ĳ������С���Ĵ�͸���
       function penetration_loss= interfering_penetration_loss(obj,interferingEnodeBids)
           pos = obj.attached_UE.pos;
           penetration_loss = reshape(obj.macroscopic_pathloss_model.get_penetration_loss_eNodeB(pos,interferingEnodeBids),[],1);
       end
       %% �õ����и���UE����UE����С���Ĵ�͸���
       function penetration_loss= ul_interfering_penetration_loss(obj,this_sector_id,interfering_ue_id)
           UE_pos_vector = obj.macroscopic_pathloss_model.UE_pos_vector;
           interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);
           for u_ = 1:size(interfering_ue_pos,1);
               penetration_loss(u_) = obj.macroscopic_pathloss_model.get_penetration_loss_eNodeB(interfering_ue_pos(u_,:),this_sector_id);
           end
       end
         %% �õ���UE������С���Ĵ������
       function pathloss_900 = macroscopic_pathloss_900(obj)
           eNodeB_id = obj.attached_UE.attached_eNodeB.eNodeB_id;
           pos = obj.attached_UE.pos;
           pathloss_900 = obj.macroscopic_pathloss_model.get_pathloss_900_eNodeB(pos,eNodeB_id);
       end
       
       %% �õ�ָ������С������UE�Ĵ������
       function pathloss_900=  interfering_macroscopic_pathloss_900(obj,interferingEnodeBids)
           pos = obj.attached_UE.pos;
           pathloss_900 = reshape(obj.macroscopic_pathloss_model.get_pathloss_900_eNodeB(pos,interferingEnodeBids),[],1);
      end
       %% �õ�����С������UE�Ĵ������
       function [pathloss is_set] = macroscopic_pathloss(obj)
           if ~obj.macroscopic_pathloss_model_is_set
               pathloss = 0;
               is_set   = false;
               return
           else

               eNodeB_id = obj.attached_UE.attached_eNodeB.eNodeB_id;
               pos = obj.attached_UE.pos;
               pathloss = obj.macroscopic_pathloss_model.get_pathloss_eNodeB(pos,eNodeB_id);
               is_set   = true;
           end
       end
       
       %% �õ�ָ������С������UE�Ĵ������
       function [pathloss is_set] = interfering_macroscopic_pathloss(obj,interferingEnodeBids)
           if ~obj.macroscopic_pathloss_model_is_set
               pathloss = zeros([length(interferingEnodeBids) 1]); % Column vector
               is_set   = false;
               return
           else
               pos = obj.attached_UE.pos;
               pathloss = reshape(obj.macroscopic_pathloss_model.get_pathloss_eNodeB(pos,interferingEnodeBids),[],1);
               is_set   = true;
           end
       end
       %% �õ����и���UE���÷����վ��·��
       function [pathloss is_set] = ul_interfering_macroscopic_pathloss(obj,this_sector_id,interfering_ue_id)
           if ~obj.macroscopic_pathloss_model_is_set
               UE_pos_vector = obj.macroscopic_pathloss_model.UE_pos_vector;
               interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);
               pathloss = zeros(1,size(interfering_ue_pos,1));
               is_set   = false;
               return
           else
               UE_pos_vector = obj.macroscopic_pathloss_model.UE_pos_vector;
               interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);
               for u_ = 1:size(interfering_ue_pos,1);
                   pathloss(u_) = obj.macroscopic_pathloss_model.get_pathloss_eNodeB(interfering_ue_pos(u_,:),this_sector_id);
               end
               is_set   = true;
           end
       end
       
       %% �õ���UE����Ӱ˥��
       function [pathloss is_set] = shadow_fading_pathloss(obj)
           if ~obj.shadow_fading_model_is_set
               pathloss = 0;
               is_set   = false;
               return
           else
               % �õ�eNodeB id ��pos
               eNodeB_id = obj.attached_UE.attached_site.id;
               pos       = obj.attached_UE.pos;
               pathloss  = obj.shadow_fading_model.get_pathloss(pos,eNodeB_id);
               is_set    = true;
           end
       end
       
       %% �õ�ָ�����Ż�վ����UE����Ӱ˥��
       function [pathloss is_set] = interfering_shadow_fading_pathloss(obj,interferingSiteIds)
           if ~obj.shadow_fading_model_is_set
               pathloss = zeros([length(interferingSiteIds) 1]); % Column vector
               is_set   = false;
               return
           else
               pos      = obj.attached_UE.pos;
               pathloss = reshape(obj.shadow_fading_model.get_pathloss(pos,interferingSiteIds),[],1);
               is_set   = true;
           end
       end
       %% �õ�����UE���÷����վ����Ӱ˥��
       function [pathloss is_set] = ul_interfering_shadow_fading_pathloss(obj,this_enb_id,interfering_ue_id)
           if ~obj.shadow_fading_model_is_set
               UE_pos_vector = obj.macroscopic_pathloss_model.UE_pos_vector;
               interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);
               pathloss = zeros(1,size(interfering_ue_pos,1));
               is_set   = false;
               return
           else
               UE_pos_vector = obj.macroscopic_pathloss_model.UE_pos_vector;
               interfering_ue_pos = UE_pos_vector(interfering_ue_id,:);
               for u_ = 1:size(interfering_ue_pos,1);
                   pathloss(u_) = obj.shadow_fading_model.get_pathloss(interfering_ue_pos(u_,:),this_enb_id);
               end
               is_set   = true;
           end
       end
       %% �õ�RB���������
       function the_RB_grid = RB_grid(obj)
           the_RB_grid = obj.attached_UE.attached_eNodeB.RB_grid;
       end
       %% ���ú�·��ģ��
       function set_macroscopic_pathloss_model(obj,macroscopic_pathloss_model)
           obj.macroscopic_pathloss_model        = macroscopic_pathloss_model;
           obj.macroscopic_pathloss_model_is_set = true;
       end
       %% ������Ӱ˥��ģ��
       function set_shadow_fading_model(obj,shadow_fading_model)
           obj.shadow_fading_model        = shadow_fading_model;
           obj.shadow_fading_model_is_set = true;
       end
   end
end 
