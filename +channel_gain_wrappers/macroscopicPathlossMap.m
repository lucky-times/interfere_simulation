classdef macroscopicPathlossMap < handle
    % ·��ͼ��
    % ��Ҫ����·��ȸ��ֲ�����Ҳ�����˵���·��ĺ���
    % ͬʱҲ�����˼�¼С���ֲ��ı����Լ����øı����ĺ��� 
    
    properties
        % ·�������������
        pathloss
        % BS����������
        BS_antenna_gain
         %UE��������
        UE_antenna_gain
        % 900ֱ�Ӽ���õ���·��
        pathloss_900
        %NLOS�ж�
        LOS_compare
        %UE�ĸ߶�
        H_height
        %UE����ת��
        H_orientation
        % ��͸���
        penetration_loss
        % ����
        distance_matrix
        micro_centre_distance_matrix
        % UEλ��
        UE_pos_vector
        % �칹��վ��
        num_hetnet_sites
        % �칹��վ������
        num_hetnet_sectors
         % ��һ��ϵͳ��վ��
        num_first_sites
        % ��һ��ϵͳ������
        num_first_sectors
        % ��һ��ϵͳUE��
        num_first_UEs
        
        % ��һ��ϵͳ��sector_assignment
        % sector_assignment��һ����ά���飬�ٸ�����sector_assignment��3,5��=11��
        % λ��3��5�е��Ǹ����ص����ķ����վidΪ11
        sector_assignment
        
        % �ڶ���ϵͳ��sector_assignment
        sector_assignment_double
        
        % ��һ��ϵͳ��sector_sizes����ʾÿ��С����������ص���
        sector_sizes
        
        % �ڶ���ϵͳ��sector_sizes
        sector_sizes_double
        
        % ��ͼ�ֱ���
        data_res
        
        % ��ͼ��XY��
        roi_x
        roi_y
        
        % �������·��ģ�͵�����
        name
        
        sector_idx_mapping  % sector_idx_mapping(s_idx)  = [b_ s_]
        site_sector_mapping % site_sector_mapping(b_,s_) = s_idx

        % �������ģ�����GUI�б��
        sector_centers
        
        % �ڶ���ϵͳ��sector_centers
        sector_centers_double
       

    end

    methods
        function print(obj)
            fprintf('macroscopicPathlossMap\n');
            fprintf('Data resolution: %d meters/pixel\n',obj.data_res);
            fprintf('ROI: x: %d,%d y:%d,%d\n',obj.roi_x(1),obj.roi_x(2),obj.roi_y(1),obj.roi_y(2));
        end
        
        %% ����һϵ�к����������Ǹ���UE���ڵ�λ�÷��ر����ڸ����е�·��pathloss��UE�߶�H_height��UE��ת�Ƕ�H_orientation��ֵ
        % �����������ƣ��ȵõ�UE�ľ���λ�ã�����ת��Ϊ���ص��λ�ã������ص��ھ����е���������Ȼ����������ھ������
        
        %������������
        function BS_antenna_gain= get_bs_antenna_gain(obj,pos,enodeB_idx)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
%             if point_outside_lower_bound || point_outside_upper_bound
%                 BS_antenna_gain = NaN;
%             else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%���ǽ�pos_to_pixel
                
               BS_antenna_gain = obj.BS_antenna_gain(pixel_coord(:,2),pixel_coord(:,1),enodeB_idx);
%             end
        end
        
 %������������
        function UE_antenna_gain= get_ue_antenna_gain(obj,pos,enodeB_idx)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                UE_antenna_gain = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%���ǽ�pos_to_pixel
                
               UE_antenna_gain = obj.UE_antenna_gain(pixel_coord(:,2),pixel_coord(:,1),enodeB_idx);
            end
        end
        
        
        % ����·��
        function pathloss = get_pathloss(obj,pos,s_,b_)
            s_idx = obj.site_sector_mapping(b_,s_);
            s_idx = s_idx(:);
            s_idx = s_idx(s_idx~=0);
            pathloss = obj.get_pathloss_eNodeB(pos,s_idx);
        end
        % ����ĳ���㵽ĳ��վ��·��
        function pathloss = get_pathloss_eNodeB(obj,pos,enodeB_idx)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                pathloss = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%���ǽ�pos_to_pixel
                
                pathloss = obj.pathloss(pixel_coord(:,2),pixel_coord(:,1),enodeB_idx);
            end
        end
        
        % ����UE�߶�
        function height = get_UE_height(obj,pos)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                height = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%���ǽ�pos_to_pixel
                
                height = obj.H_height(pixel_coord(:,2),pixel_coord(:,1));
            end
        end
        
        % ����UE��ת�Ƕ�
         function orientation = get_UE_orientation(obj,pos)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                orientation = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%���ǽ�pos_to_pixel
                
                orientation = obj.H_orientation(pixel_coord(:,2),pixel_coord(:,1));
            end
         end
        
        % ����û�����������·��
        function pathloss_900 = get_pathloss_900_eNodeB(obj,pos,enodeB_idx)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                pathloss_900 = NaN;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%���ǽ�pos_to_pixel
                
                pathloss_900 = obj.pathloss_900(pixel_coord(:,2),pixel_coord(:,1),enodeB_idx);
            end
        end
        
        % ���ش�͸���
        function penetration_loss = get_penetration_loss_eNodeB(obj,pos,enodeB_idx)
            x_ = pos(1);
            y_ = pos(2);
            point_outside_lower_bound = sum([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]);
            point_outside_upper_bound = sum([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]);
            if point_outside_lower_bound || point_outside_upper_bound
                penetration_loss = 0;
            else
                pixel_coord(:,1) = floor((x_-obj.roi_x(1))/obj.data_res)+1;
                pixel_coord(:,2) = floor((y_-obj.roi_y(1))/obj.data_res)+1;%���ǽ�pos_to_pixel
                
                penetration_loss = obj.penetration_loss(pixel_coord(:,2),pixel_coord(:,1),enodeB_idx);
            end
        end
        
        % Plots the pathloss of a given eNodeB sector
        function plot_pathloss(obj,b_,s_)
            figure;
            s_idx = obj.site_sector_mapping(b_,s_);
            imagesc(obj.pathloss(:,:,s_idx));
        end
        
        % Range of positions in which there are valid pathloss values
        function [x_range y_range] = valid_range(obj)
            x_range = [ obj.roi_x(1) obj.roi_x(2) ];
            y_range = [ obj.roi_y(1) obj.roi_y(2) ];
        end
        
        % Returns the coordinate origin for this pathloss map
        function pos = coordinate_origin(obj)
            pos = [ obj.roi_x(1) obj.roi_y(1) ];
        end
        
        % ��LTEϵͳ����MCL
        function apply_MCL(obj,sectors,minimum_coupling_loss)
            obj.pathloss(:,:,sectors) = max(obj.pathloss(:,:,sectors),minimum_coupling_loss);
        end
        
        % ֪��UE��λ�ã����ظ�UE��������������վ
        function [ b_ s_ s_idx] = cell_assignment(obj,SYS_config,UE_id,pos,num_first_UEs,varargin)
            if length(varargin)<1
                x_ = pos(1);
                y_ = pos(2);
            else
                x_ = pos(1);
                y_ = varargin{1};
            end
            if max([x_ y_] < [obj.roi_x(1),obj.roi_y(1)]) >= 1
                eNodeB_id = NaN;
                sector_num = NaN;
            elseif max([x_ y_] > [obj.roi_x(2),obj.roi_y(2)]) >= 1
                eNodeB_id = NaN;
                sector_num = NaN;
            else
                pixel_coord = NR_common_pos_to_pixel([x_ y_],obj.coordinate_origin,obj.data_res);
                if SYS_config.isDouble
                    if UE_id<=num_first_UEs
                        s_idx = obj.sector_assignment(pixel_coord(2),pixel_coord(1));
                    else
                        s_idx = obj.sector_assignment_double(pixel_coord(2),pixel_coord(1));
                    end
                else
                    s_idx = obj.sector_assignment(pixel_coord(2),pixel_coord(1));%�ҵ���UE����λ�ö�Ӧ��������Ȼ�����ҵ���Ӧ��enb
                end
                
                b_s_        = obj.sector_idx_mapping(s_idx,:);
                b_          = b_s_(1);
                s_          = b_s_(2);
            end
        end
        
        % ���ػ�վ����������
        function [ eNodeBs sectors_per_eNodeB ] = size(obj)
            eNodeBs            = size(obj.pathloss,4);
            sectors_per_eNodeB = size(obj.pathloss,3);
        end
        
        % ��UE�˶��Ƴ��˵�ͼ�����¸�������һ��λ��
        function position = random_position(obj)
            [row,col] = find(obj.sector_assignment > 0);
            index = NR_randperm(length(col),1);
            position_pix = [col(index),row(index)];
            position = NR_common_pixel_to_pos(position_pix,obj.coordinate_origin,obj.data_res);
        end
        
        % Deletes all information except for the cell assignment maps
        function delete_everything_except_cell_assignments(obj)
            obj.pathloss            = [];
            obj.name                = [];
            obj.sector_idx_mapping  = [];
            obj.site_sector_mapping = [];
            obj.BS_antenna_gain     = [];
            obj.UE_antenna_gain     = [];
            obj.pathloss_900        = [];
            obj.LOS_compare         = [];
            obj.H_height            = [];
            obj.H_orientation       = [];
            obj.penetration_loss    = [];
            obj.distance_matrix     = [];
            obj.micro_centre_distance_matrix  = [];
            obj.UE_pos_vector       = [];
        end
    end
end