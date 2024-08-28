function [ eNodeB_sites,num_hetnet_sites] = NR_init_create_eNodeBs(SYS_config)
% ������վ��λ��
% ���������
% SYS_config����������
%
%���������
%eNodeB_sites����վʵ����Ϣ
% num_hetnet_sites���칹�����У��ڶ���ϵͳ��վ��

num_hetnet_sites = 0;% ���칹�����⣬��ֵ��Ϊ0

%��Բ�ͬ����������Ӧ����
switch SYS_config.scene_type
    case {'UMA','RMa'}
%         n_rings = SYS_config.nr_eNodeB_rings;% UMa�л�վΧ������Ȧ����0����ֻ��һ����վ
%         shift_mode = SYS_config.shift_mode;% ƫ��
%         
        ISD = SYS_config.ISD;% ��վ�����
%         [tmp_gridx,tmp_gridy] = meshgrid(-n_rings:n_rings,...
%             (-n_rings:n_rings)*sin(pi/3));
%         % �����㣬����������������ʽ
%         % * * * * *
%         % * * * * *
%         % * * * * *
%         % * * * * *
%         % * * * * *
%         if mod(n_rings,2) == 0
%             tmp_shift_idx = 2:2:2*n_rings+1; %ת��ż����
%             % ����
%             % * * * * *
%             %  * * * * *
%             % * * * * *
%             %  * * * * *
%             % * * * * *
%         else
%             tmp_shift_idx = 1:2:2*n_rings+1; %ת��������
%             % ����
%             %  * * * * *
%             % * * * * *
%             %  * * * * *
%             % * * * * *
%             %  * * * * *
%         end
%         
%         tmp_gridx(tmp_shift_idx,:) = tmp_gridx(tmp_shift_idx,:) + 0.5;%ת��
%         
%         rot = @(w_) [cos(w_),-sin(w_);sin(w_),cos(w_)]; %rot��һ��������w_������������������������ǲ���һ������ ˳ʱ����תw_
%         for i_ = 1:7
%             tmp_hex(i_,:) = ((n_rings+0.5)*rot(pi/3)^(i_-1)*[1;0]).'; % ��[1;0]��ʾֻȡ��һ��
%         end
%         
%         tmp_valid_positions = inpolygon(tmp_gridx,tmp_gridy,tmp_hex(:,1),tmp_hex(:,2));% tmp_hex���������ζ������꣬inpolygon�жϵ��Ƿ�����������
%         tmp_x = tmp_gridx(tmp_valid_positions);
%         tmp_y = tmp_gridy(tmp_valid_positions);

        %���ݸ�����excel�����BSվ��ľ�γ�ȣ�ȷ������BS��վ��λ��
        data2 = readtable('7��4.9Gվ��.xlsx');
        info = data2(:, [9, 10]);
        info = table2array(info);
        pos = [];
        for i = 1:length(info)
            if(ismember(info(i, 1), pos))
                continue;
            else
                pos = [pos;info(i, :)];
            end
        end
%         figure;
%         hold on;
        for b_ = 1:length(pos)
            eNodeB_sites(b_)           = network_elements.eNodeB_site;
            eNodeB_sites(b_).id        = b_;
%             eNodeB_sites(b_).pos       = [tmp_x(b_)*ISD tmp_y(b_)*ISD];
            eNodeB_sites(b_).pos       = latlon_to_xy(pos(b_, :));
            pos_eNodeB(b_, :) = eNodeB_sites(b_).pos;
            SYS_config.eNodeB_pos(b_, :) = eNodeB_sites(b_).pos;
%             scatter(eNodeB_sites(b_).pos(1), eNodeB_sites(b_).pos(2));
%             text(eNodeB_sites(b_).pos(1), eNodeB_sites(b_).pos(2), num2str(b_), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
            eNodeB_sites(b_).site_type = 'macro';
            % �����������ģ����ڻ�����
            sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+ISD/3 eNodeB_sites(b_).pos(2)];% ���㷽���Ǹ��������ε��������� ���ҷ� ISD/3
            sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];  %���� ISD/3 ��ʱ����ת30��
            sector_centre(3).pos = [eNodeB_sites(b_).pos(1)-(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];%���ҷ�˳ʱ����ת60�� ISD/3 ˳ʱ����ת30��
            eNodeB_sites(b_).sector_centre = sector_centre;
        end
%         scatter(pos_eNodeB(:, 1), pos_eNodeB(:, 2), "filled", 'r');
        % ��һ��ϵͳ��wraparound��վ
        if SYS_config.isWraparound 
            current_enb = b_;
            n_rings = 4;%������Ȧ
            [tmp_gridx,tmp_gridy] = meshgrid(-n_rings:n_rings,...
                (-n_rings:n_rings)*sin(pi/3)); 
            if mod(n_rings,2) == 0
                tmp_shift_idx = 2:2:2*n_rings+1;
            else
                tmp_shift_idx = 1:2:2*n_rings+1; 
            end
            
            tmp_gridx(tmp_shift_idx,:) = tmp_gridx(tmp_shift_idx,:) + 0.5; 
            
            rot = @(w_) [cos(w_),-sin(w_);sin(w_),cos(w_)]; 
            for i_ = 1:7
               
                tmp_hex(i_,:) = ((n_rings+0.5)*rot(pi/3)^(i_-1)*[1;0]).'; 
            end
            
            tmp_valid_positions = inpolygon(tmp_gridx,tmp_gridy,tmp_hex(:,1),tmp_hex(:,2));
            tmp_x_1 = tmp_gridx(tmp_valid_positions);
            tmp_y_1 = tmp_gridy(tmp_valid_positions);
            % ��������������ֱ�������֣�ֻ�����ַ���������ʹ��Ԫ�����齫����ת�����ַ���
            tmp = cell(1,length(tmp_x));
            tmp1 = cell(1,length(tmp_x_1));            
            for i = 1:length(tmp_x)
                tmp{i} = [num2str(tmp_x(i)),num2str(tmp_y(i))];
            end
            for i = 1:length(tmp_x_1)
                tmp1{i} = [num2str(tmp_x_1(i)),num2str(tmp_y_1(i))];
            end            
            [~,~,ib] = intersect(tmp,tmp1);
            
            % �½���һ��5Ȧ�����ˣ�Ȼ����ԭ����������������5Ȧ�����м�3Ȧȥ����
            % �����������˺����������м�3ȦidС��������Ȧid�������
            tmp_x_1(ib) = [];
            tmp_y_1(ib) = [];
            
            for b_ = current_enb+1:current_enb+length(tmp_x_1);
                eNodeB_sites(b_)           = network_elements.eNodeB_site;
                eNodeB_sites(b_).id        = b_;
                eNodeB_sites(b_).pos       = [tmp_x_1(b_-current_enb)*ISD tmp_y_1(b_-current_enb)*ISD];
                eNodeB_sites(b_).site_type = 'macro';
                % ������������
                sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];
                sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-ISD/3 eNodeB_sites(b_).pos(2)];
                sector_centre(3).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];
                eNodeB_sites(b_).sector_centre = sector_centre;
            end
        end
        
        % ���ݵ�ϵͳ/˫ϵͳ�ֱ���Թ�ַ/ƫ�ò�������
        if SYS_config.isDouble %˫ϵͳ
            current_enb1 = b_;
            switch shift_mode
                case 0 %��ַ
                    for b_ = current_enb1+1:current_enb1+length(tmp_x)
                        eNodeB_sites(b_)           = network_elements.eNodeB_site;
                        eNodeB_sites(b_).id        = b_;
                        eNodeB_sites(b_).pos       = [tmp_x(b_-current_enb1)*ISD tmp_y(b_-current_enb1)*ISD];
                        eNodeB_sites(b_).site_type = 'macro';
                        % ������������
                        sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];
                        sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-ISD/3 eNodeB_sites(b_).pos(2)];
                        sector_centre(3).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];
                        eNodeB_sites(b_).sector_centre = sector_centre;
                    end
                case 1 %ƫ��
                    for b_ = current_enb1+1:current_enb1+length(tmp_x)
                        eNodeB_sites(b_)           = network_elements.eNodeB_site;
                        eNodeB_sites(b_).id        = b_;
                        eNodeB_sites(b_).pos       = [tmp_x(b_-current_enb1)*ISD-(ISD/2) tmp_y(b_-current_enb1)*ISD-(ISD*(3^0.5)/6)];
                        eNodeB_sites(b_).site_type = 'macro';
                        % ������������
                        sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];
                        sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-ISD/3 eNodeB_sites(b_).pos(2)];
                        sector_centre(3).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];
                        eNodeB_sites(b_).sector_centre = sector_centre;
                    end
                case 2 % �������õ�ƫ��
                    for b_ = current_enb1+1:current_enb1+length(tmp_x)
                        eNodeB_sites(b_)           = network_elements.eNodeB_site;
                        eNodeB_sites(b_).id        = b_;
                        eNodeB_sites(b_).pos       = [tmp_x(b_-current_enb1)*ISD+SYS_config.ISD_between_systems*cos(SYS_config.angle_between_systems/180*pi) tmp_y(b_-current_enb1)*ISD+SYS_config.ISD_between_systems*sin(SYS_config.angle_between_systems/180*pi)];
                        eNodeB_sites(b_).site_type = 'macro';
                        % ������������
                        sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];
                        sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-ISD/3 eNodeB_sites(b_).pos(2)];
                        sector_centre(3).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];
                        eNodeB_sites(b_).sector_centre = sector_centre;
                    end
            end
            
            % �ڶ���ϵͳ��wraparound
            if SYS_config.isWraparound
                current_enb2 = b_;
                switch shift_mode
                    case 0
                        for b_ = current_enb2+1:current_enb2+length(tmp_x_1);
                            eNodeB_sites(b_)           = network_elements.eNodeB_site;
                            eNodeB_sites(b_).id        = b_;
                            eNodeB_sites(b_).pos       = [tmp_x_1(b_-current_enb2)*ISD tmp_y_1(b_-current_enb2)*ISD];
                            eNodeB_sites(b_).site_type = 'macro';
                            % ������������
                            sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];
                            sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-ISD/3 eNodeB_sites(b_).pos(2)];
                            sector_centre(3).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];
                            eNodeB_sites(b_).sector_centre = sector_centre;
                        end
                    case 1
                        for b_ = current_enb2+1:current_enb2+length(tmp_x_1);
                            eNodeB_sites(b_)           = network_elements.eNodeB_site;
                            eNodeB_sites(b_).id        = b_;
                            eNodeB_sites(b_).pos       = [tmp_x_1(b_-current_enb2)*ISD-(ISD/2) tmp_y_1(b_-current_enb2)*ISD-(ISD*(3^0.5)/6)];
                            eNodeB_sites(b_).site_type = 'macro';
                            % ������������
                            sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];
                            sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-ISD/3 eNodeB_sites(b_).pos(2)];
                            sector_centre(3).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];
                            eNodeB_sites(b_).sector_centre = sector_centre;
                        end
                    case 2
                        for b_ = current_enb2+1:current_enb2+length(tmp_x_1);
                            eNodeB_sites(b_)           = network_elements.eNodeB_site;
                            eNodeB_sites(b_).id        = b_;
                            eNodeB_sites(b_).pos       = [tmp_x_1(b_-current_enb2)*ISD+SYS_config.ISD_between_systems*cos(SYS_config.angle_between_systems/180*pi) tmp_y_1(b_-current_enb2)*ISD+SYS_config.ISD_between_systems*sin(SYS_config.angle_between_systems/180*pi)];
                            eNodeB_sites(b_).site_type = 'macro';
                            % ������������
                            sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];
                            sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-ISD/3 eNodeB_sites(b_).pos(2)];
                            sector_centre(3).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];
                            eNodeB_sites(b_).sector_centre = sector_centre;
                        end
                end
            end
            
        end
    case {'UMI','UMa_to_UMi'} %UMI��������UMA,UMI�칹����
        if ~SYS_config.isManhattan %����Ϊmanhattan����
            n_rings = SYS_config.nr_eNodeB_rings;
            ISD = SYS_config.ISD;
            UMi_r = SYS_config.UMi_r;
            % ȷ������ѵ�λ��
            % �����UMa�Ļ���һ��
            [tmp_gridx,tmp_gridy] = meshgrid(-n_rings:n_rings,...
                (-n_rings:n_rings)*sin(pi/3)); 
            if mod(n_rings,2) == 0
                tmp_shift_idx = 2:2:2*n_rings+1; 
            else
                tmp_shift_idx = 1:2:2*n_rings+1; 
            end
            
            tmp_gridx(tmp_shift_idx,:) = tmp_gridx(tmp_shift_idx,:) + 0.5;
            
            rot = @(w_) [cos(w_),-sin(w_);sin(w_),cos(w_)]; 
            for i_ = 1:7
           
                tmp_hex(i_,:) = ((n_rings+0.5)*rot(pi/3)^(i_-1)*[1;0]).'; 
            end
            
            tmp_valid_positions = inpolygon(tmp_gridx,tmp_gridy,tmp_hex(:,1),tmp_hex(:,2));
            tmp_x = tmp_gridx(tmp_valid_positions);
            tmp_y = tmp_gridy(tmp_valid_positions);
            
            for vb_ = 1:length(tmp_x)
                veNodeBs(vb_).pos       = [tmp_x(vb_)*ISD tmp_y(vb_)*ISD];
            end
            
            % ����ĺ��վ����������վ����Ϊ���ݳ�����ȷ���Ƿ���Ҫ����MBS
            % ������random dropֻ��Ҫ������MBS��ȷ��΢��վλ��
            vb_ = 1;
            for s_ = 1:3:3*length(tmp_x)
                sector_centre(s_).pos = [veNodeBs(vb_).pos(1)+(ISD/3)/2 veNodeBs(vb_).pos(2)+(ISD/3)/2*(3^0.5)];
                sector_centre(s_+1).pos = [veNodeBs(vb_).pos(1)-ISD/3 veNodeBs(vb_).pos(2)];
                sector_centre(s_+2).pos = [veNodeBs(vb_).pos(1)+(ISD/3)/2 veNodeBs(vb_).pos(2)-(ISD/3)/2*(3^0.5)];
                vb_ = vb_+1;
            end
            
            b_ = 1;
            t = 0:360;
            theta = 2*pi/360 * t;
            % ȷ��΢��վ�Ĳ���Χ
            for s_ = 1:3*length(tmp_x)
                micro_centre = [sector_centre(s_).pos(1)-(ISD/3)/2/2 sector_centre(s_).pos(2)+(ISD/3)/2/2*(3^0.5)];% ��С���������΢С������
                tmp = micro_centre(1) + UMi_r * cos(theta);% ������������õ�΢С����վλ��
                y = micro_centre(2) + UMi_r * sin(theta);
                index = NR_randperm(361,1);
                eNodeB_sites(b_)           = network_elements.eNodeB_site;
                eNodeB_sites(b_).id        = b_;
                eNodeB_sites(b_).pos       = [tmp(index) y(index)];
                eNodeB_sites(b_).parent_centre_pos = micro_centre;
                eNodeB_sites(b_).site_type = 'micro';
                eNodeB_sites(b_).sector_centre = sector_centre(s_);
                b_ = b_+1;
                
                micro_centre = [sector_centre(s_).pos(1)-(ISD/3)/2/2 sector_centre(s_).pos(2)-(ISD/3)/2/2*(3^0.5)];
                tmp = micro_centre(1) + UMi_r * cos(theta);
                y = micro_centre(2) + UMi_r * sin(theta);
                index = NR_randperm(361,1);
                eNodeB_sites(b_)           = network_elements.eNodeB_site;
                eNodeB_sites(b_).id        = b_;
                eNodeB_sites(b_).pos       = [tmp(index) y(index)];
                eNodeB_sites(b_).parent_centre_pos = micro_centre;
                eNodeB_sites(b_).site_type = 'micro';
                eNodeB_sites(b_).sector_centre = sector_centre(s_);
                b_ = b_+1;
                
                micro_centre = [sector_centre(s_).pos(1)+(ISD/3)/2 sector_centre(s_).pos(2)];
                tmp = micro_centre(1) + UMi_r * cos(theta);
                y = micro_centre(2) + UMi_r * sin(theta);
                index = NR_randperm(361,1);
                eNodeB_sites(b_)           = network_elements.eNodeB_site;
                eNodeB_sites(b_).id        = b_;
                eNodeB_sites(b_).pos       = [tmp(index) y(index)];
                eNodeB_sites(b_).parent_centre_pos = micro_centre;
                eNodeB_sites(b_).site_type = 'micro';
                eNodeB_sites(b_).sector_centre = sector_centre(s_);
                b_ = b_+1;
            end
            
            % ˫ϵͳ
            if SYS_config.isDouble
                switch SYS_config.scene_type
                    case 'UMI'
                        for s_ = 1:3*length(tmp_x)
                            micro_centre = [sector_centre(s_).pos(1)-(ISD/3)/2/2 sector_centre(s_).pos(2)+(ISD/3)/2/2*(3^0.5)];
                            tmp = micro_centre(1) + UMi_r * cos(theta);
                            y = micro_centre(2) + UMi_r * sin(theta);
                            index = NR_randperm(361,1);
                            eNodeB_sites(b_)           = network_elements.eNodeB_site;
                            eNodeB_sites(b_).id        = b_;
                            eNodeB_sites(b_).pos       = [tmp(index) y(index)];
                            eNodeB_sites(b_).parent_centre_pos = micro_centre;
                            eNodeB_sites(b_).site_type = 'micro';
                            eNodeB_sites(b_).sector_centre = sector_centre(s_);
                            b_ = b_+1;
                            
                            micro_centre = [sector_centre(s_).pos(1)-(ISD/3)/2/2 sector_centre(s_).pos(2)-(ISD/3)/2/2*(3^0.5)];
                            tmp = micro_centre(1) + UMi_r * cos(theta);
                            y = micro_centre(2) + UMi_r * sin(theta);
                            index = NR_randperm(361,1);
                            eNodeB_sites(b_)           = network_elements.eNodeB_site;
                            eNodeB_sites(b_).id        = b_;
                            eNodeB_sites(b_).pos       = [tmp(index) y(index)];
                            eNodeB_sites(b_).parent_centre_pos = micro_centre;
                            eNodeB_sites(b_).site_type = 'micro';
                            eNodeB_sites(b_).sector_centre = sector_centre(s_);
                            b_ = b_+1;
                            
                            micro_centre = [sector_centre(s_).pos(1)+(ISD/3)/2 sector_centre(s_).pos(2)];
                            tmp = micro_centre(1) + UMi_r * cos(theta);
                            y = micro_centre(2) + UMi_r * sin(theta);
                            index = NR_randperm(361,1);
                            eNodeB_sites(b_)           = network_elements.eNodeB_site;
                            eNodeB_sites(b_).id        = b_;
                            eNodeB_sites(b_).pos       = [tmp(index) y(index)];
                            eNodeB_sites(b_).parent_centre_pos = micro_centre;
                            eNodeB_sites(b_).site_type = 'micro';
                            eNodeB_sites(b_).sector_centre = sector_centre(s_);
                            b_ = b_+1;
                        end
                    case 'UMa_to_UMi'
                        current_eNodeBs = length(eNodeB_sites);
                        num_hetnet_sites = length(veNodeBs);
                        for b_ = current_eNodeBs+1:current_eNodeBs+length(veNodeBs)
                            eNodeB_sites(b_)           = network_elements.eNodeB_site;
                            eNodeB_sites(b_).id        = b_;
                            eNodeB_sites(b_).pos       = veNodeBs(b_-current_eNodeBs).pos;% ��������õ�������վλ�ø��������ĺ��վ
                            eNodeB_sites(b_).site_type = 'macro';
                            % ������������
                            sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];
                            sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-ISD/3 eNodeB_sites(b_).pos(2)];
                            sector_centre(3).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];
                            eNodeB_sites(b_).sector_centre = sector_centre;
                        end
                end
            end
        else %�����manhattan����
            % �����ٵĻ�վλ�ö��ǹ̶���
            % pos_x�����л�վ��x��
            pos_x = [217.5 577.5 937.5 82.5 442.5 802.5 307.5 667.5 1027.5 172.5 532.5 892.5 37.5 397.5 757.5 262.5 622.5 982.5 127.5 487.5 847.5 352.5 712.5 1072.5];
            index = 1;
            row = 1;
            for i=1:3
                for b_ = (i-1)*24+1:i*24;
                    eNodeB_sites(b_)           = network_elements.eNodeB_site;
                    eNodeB_sites(b_).id        = b_;
                    eNodeB_sites(b_).pos       = [pos_x(index) 7.5+(row-1)*45];% index��row��Ϊ����ȷ����վ��x��y��
                    switch SYS_config.scene_type
                        case 'UMa_to_UMi'
                            eNodeB_sites(b_).pos = [eNodeB_sites(b_).pos(1)-540,eNodeB_sites(b_).pos(2)-540];
                    end
                    index = index+1;
                    if mod(b_,3) == 0
                        row = row+1;
                    end
                    eNodeB_sites(b_).site_type = 'micro';
                end
                index = 1;
            end
            if SYS_config.isDouble
                switch SYS_config.scene_type
                    case 'UMI'
                        pos_x = [37.5 397.5 757.5 262.5 622.5 982.5 127.5 487.5 847.5 352.5 712.5 1072.5 217.5 577.5 937.5 82.5 442.5 802.5 307.5 667.5 1027.5 172.5 532.5 892.5];
                        index = 1;
                        row = 1;
                        for i=4:6
                            for b_ = (i-1)*24+1:i*24;
                                eNodeB_sites(b_)           = network_elements.eNodeB_site;
                                eNodeB_sites(b_).id        = b_;
                                eNodeB_sites(b_).pos       = [pos_x(index) 7.5+(row-1)*45];
                                index = index+1;
                                if mod(b_,3) == 0
                                    row = row+1;
                                end
                                eNodeB_sites(b_).site_type = 'micro';
                            end
                            index = 1;
                        end
                    case 'UMa_to_UMi'
                        current_eNodeBs = length(eNodeB_sites);
                        n_rings = SYS_config.nr_eNodeB_rings;
                        ISD = SYS_config.ISD;
                        % ȷ������ѵ�λ��
                        % �����UMa�Ļ���һ��
                        [tmp_gridx,tmp_gridy] = meshgrid(-n_rings:n_rings,...
                            (-n_rings:n_rings)*sin(pi/3));
                        if mod(n_rings,2) == 0
                            tmp_shift_idx = 2:2:2*n_rings+1;
                        else
                            tmp_shift_idx = 1:2:2*n_rings+1;
                        end
                        
                        tmp_gridx(tmp_shift_idx,:) = tmp_gridx(tmp_shift_idx,:) + 0.5;
                        
                        rot = @(w_) [cos(w_),-sin(w_);sin(w_),cos(w_)];
                        for i_ = 1:7
                            tmp_hex(i_,:) = ((n_rings+0.5)*rot(pi/3)^(i_-1)*[1;0]).';
                        end
                        
                        tmp_valid_positions = inpolygon(tmp_gridx,tmp_gridy,tmp_hex(:,1),tmp_hex(:,2));
                        tmp_x = tmp_gridx(tmp_valid_positions);
                        tmp_y = tmp_gridy(tmp_valid_positions);
                        for b_ = current_eNodeBs+1:current_eNodeBs+length(tmp_x);
                            eNodeB_sites(b_)           = network_elements.eNodeB_site;
                            eNodeB_sites(b_).id        = b_;
                            eNodeB_sites(b_).pos       = [tmp_x(b_-current_eNodeBs)*ISD tmp_y(b_-current_eNodeBs)*ISD];
                            if SYS_config.shift_mode == 2
                                eNodeB_sites(b_).pos = [eNodeB_sites(b_).pos(1)+SYS_config.ISD_between_systems*cos(SYS_config.angle_between_systems/180*pi) eNodeB_sites(b_).pos(2)+SYS_config.ISD_between_systems*sin(SYS_config.angle_between_systems/180*pi)];
                            end
                            eNodeB_sites(b_).site_type = 'macro';
                            % ������������
                            sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];
                            sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-ISD/3 eNodeB_sites(b_).pos(2)];
                            sector_centre(3).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];
                            eNodeB_sites(b_).sector_centre = sector_centre;
                        end
                        current_eNodeBs2 = length(eNodeB_sites);
                        num_hetnet_sites = current_eNodeBs2 - current_eNodeBs;
                end
            end
        end
    case 'InH' %�����InH����
        if SYS_config.isFemto
            % femto��һ��120*40��С
            row1 = 10;% �����������ڵ��У���վλ�ò�������ȷ��
            row2 = 30;
        else
            % open office ��һ��120*50��С
            row1 = 15;% ��վ���ڵ���
            row2 = 35;
        end
        for b_ = 1:12;
            eNodeB_sites(b_)           = network_elements.eNodeB_site;
            eNodeB_sites(b_).id        = b_;
            if b_<=6
                eNodeB_sites(b_).pos       = [10+(b_-1)*20 row1];% �̶�y�ᣬ�ı�x��λ��
            elseif b_>=7 && b_<=12
                eNodeB_sites(b_).pos       = [10+(b_-7)*20 row2];
            end
            eNodeB_sites(b_).site_type = 'indoor';
        end
        if SYS_config.isDouble
            for b_ = 13:24;
                eNodeB_sites(b_)           = network_elements.eNodeB_site;
                eNodeB_sites(b_).id        = b_;
                if b_<=18
                    eNodeB_sites(b_).pos       = [10+(b_-13)*20 row1];
                elseif b_>=19 && b_<=24
                    eNodeB_sites(b_).pos       = [10+(b_-19)*20 row2];
                end
                eNodeB_sites(b_).site_type = 'indoor';
            end
        end
    case 'InH2' %�����InH����
        
            % open office ��һ��120*50��С
            row1 = 15;% ��վ���ڵ���
            row2 = 35;

        for b_ = 1:6
            eNodeB_sites(b_)           = network_elements.eNodeB_site;
            eNodeB_sites(b_).id        = b_;
            if b_<=3
                eNodeB_sites(b_).pos       = [10+(b_-1)*40 row1];% �̶�y�ᣬ�ı�x��λ��
            elseif b_>=3 && b_<=6
                eNodeB_sites(b_).pos       = [30+(b_-3)*40 row2];
            end
            eNodeB_sites(b_).site_type = 'indoor';
        end
        if SYS_config.isDouble
            for b_ = 7:12
                eNodeB_sites(b_)           = network_elements.eNodeB_site;
                eNodeB_sites(b_).id        = b_;
                if b_<=9
                    eNodeB_sites(b_).pos       = [30+(b_-7)*40 row1];
                elseif b_>=10 && b_<=12
                    eNodeB_sites(b_).pos       = [10+(b_-10)*40 row2];
                end
                eNodeB_sites(b_).site_type = 'indoor';
            end
        end
    case 'InH3' %�����InH����
        
        % open office ��һ��120*50��С
        row1 = 15;% ��վ���ڵ���
        row2 = 35;
        
        for b_ = 1:6
            eNodeB_sites(b_)           = network_elements.eNodeB_site;
            eNodeB_sites(b_).id        = b_;
            if b_<=3
                eNodeB_sites(b_).pos       = [10+(b_-1)*40 row1];% �̶�y�ᣬ�ı�x��λ��
            elseif b_>=3 && b_<=6
                eNodeB_sites(b_).pos       = [30+(b_-3)*40 row2];
            end
            eNodeB_sites(b_).site_type = 'indoor';
        end
        if SYS_config.isDouble
            for b_ = 7:12
                eNodeB_sites(b_)           = network_elements.eNodeB_site;
                eNodeB_sites(b_).id        = b_;
                if b_<=9
                    eNodeB_sites(b_).pos       = [10+(b_-7)*40 row1+SYS_config.InH3_d+50];
                elseif b_>=10 && b_<=12
                    eNodeB_sites(b_).pos       = [30+(b_-10)*40 row2+SYS_config.InH3_d+50];
                end
                eNodeB_sites(b_).site_type = 'indoor';
            end
        end
    case {'UMa_to_InH','UMi_to_InH'}% Uma����Inh��Umi����Inh
        O = [0,0];% ��ͼ��ԭ��
        ISD = SYS_config.ISD;
        if ~SYS_config.isFemto
            O_InH = [O(1)-60,O(2)-25];% InH������ԭ��
            row1 = 15;
            row2 = 35;
        else
            O_InH = [O(1)-60,O(2)-20];% InH������ԭ��
            row1 = 10;
            row2 = 30;
        end
        for b_ = 1:12;
            eNodeB_sites(b_)           = network_elements.eNodeB_site;
            eNodeB_sites(b_).id        = b_;
            if b_<=6
                eNodeB_sites(b_).pos       = O_InH + [10+(b_-1)*20 row1];
            elseif b_>=7 && b_<=12
                eNodeB_sites(b_).pos       = O_InH + [10+(b_-7)*20 row2];
            end
            eNodeB_sites(b_).site_type = 'indoor';
            eNodeB_sites(b_).O_InH = O_InH;% ��������������Ժ����
        end
        
        % ������վ��λ��
        tmp(1).pos = [O(1),O(2)+(sqrt(3)/2-sqrt(3)/6)*ISD];
        tmp(2).pos = [O(1)-ISD/2,O(2)-sqrt(3)/6*ISD];
        tmp(3).pos = [O(1)+ISD/2,O(2)-sqrt(3)/6*ISD];
        % ����ϵͳ��ƫ��
        if SYS_config.shift_mode == 2
            tmp(1).pos = [tmp(1).pos(1)+SYS_config.ISD_between_systems*cos(SYS_config.angle_between_systems/180*pi) tmp(1).pos(2)+SYS_config.ISD_between_systems*sin(SYS_config.angle_between_systems/180*pi)];
            tmp(2).pos = [tmp(2).pos(1)+SYS_config.ISD_between_systems*cos(SYS_config.angle_between_systems/180*pi) tmp(2).pos(2)+SYS_config.ISD_between_systems*sin(SYS_config.angle_between_systems/180*pi)];
            tmp(3).pos = [tmp(3).pos(1)+SYS_config.ISD_between_systems*cos(SYS_config.angle_between_systems/180*pi) tmp(3).pos(2)+SYS_config.ISD_between_systems*sin(SYS_config.angle_between_systems/180*pi)];
        end
        
        current_eNodeBs = length(eNodeB_sites);% ��һ��ϵͳ�Ļ�վ��
        % ���վ��΢��վ�Ĺ�����UMa��UMi�����Ĳ��
        switch SYS_config.scene_type
            case 'UMa_to_InH'
                for b_ = current_eNodeBs+1:current_eNodeBs+3
                    eNodeB_sites(b_)           = network_elements.eNodeB_site;
                    eNodeB_sites(b_).id        = b_;
                    eNodeB_sites(b_).pos       = tmp(b_-current_eNodeBs).pos;
                    eNodeB_sites(b_).site_type = 'macro';
                    eNodeB_sites(b_).O_InH = O_InH;
                    % ������������
                    sector_centre(1).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)+(ISD/3)/2*(3^0.5)];
                    sector_centre(2).pos = [eNodeB_sites(b_).pos(1)-ISD/3 eNodeB_sites(b_).pos(2)];
                    sector_centre(3).pos = [eNodeB_sites(b_).pos(1)+(ISD/3)/2 eNodeB_sites(b_).pos(2)-(ISD/3)/2*(3^0.5)];
                    eNodeB_sites(b_).sector_centre = sector_centre;
                end
            case 'UMi_to_InH'
                UMi_r = SYS_config.UMi_r;
                vb_ = 1;
                for s_ = 1:3:3*length(tmp)
                    sector_centre(s_).pos = [tmp(vb_).pos(1)+(ISD/3)/2 tmp(vb_).pos(2)+(ISD/3)/2*(3^0.5)];
                    sector_centre(s_+1).pos = [tmp(vb_).pos(1)-ISD/3 tmp(vb_).pos(2)];
                    sector_centre(s_+2).pos = [tmp(vb_).pos(1)+(ISD/3)/2 tmp(vb_).pos(2)-(ISD/3)/2*(3^0.5)];
                    vb_ = vb_+1;
                end
                b_ = current_eNodeBs + 1;
                t = 0:360;
                theta = 2*pi/360 * t;
                % ȷ��΢��վ�Ĳ���Χ
                for s_ = 1:3*length(tmp)
                    micro_centre = [sector_centre(s_).pos(1)-(ISD/3)/2/2 sector_centre(s_).pos(2)+(ISD/3)/2/2*(3^0.5)];
                    tmp = micro_centre(1) + UMi_r * cos(theta);
                    y = micro_centre(2) + UMi_r * sin(theta);
                    index = randperm(361,1);
                    eNodeB_sites(b_)           = network_elements.eNodeB_site;
                    eNodeB_sites(b_).id        = b_;
                    eNodeB_sites(b_).pos       = [tmp(index) y(index)];
                    eNodeB_sites(b_).parent_centre_pos = micro_centre;
                    eNodeB_sites(b_).site_type = 'micro';
                    eNodeB_sites(b_).sector_centre = sector_centre(s_);
                    eNodeB_sites(b_).O_InH = O_InH;
                    b_ = b_+1;
                    
                    micro_centre = [sector_centre(s_).pos(1)-(ISD/3)/2/2 sector_centre(s_).pos(2)-(ISD/3)/2/2*(3^0.5)];
                    tmp = micro_centre(1) + UMi_r * cos(theta);
                    y = micro_centre(2) + UMi_r * sin(theta);
                    index = randperm(361,1);
                    eNodeB_sites(b_)           = network_elements.eNodeB_site;
                    eNodeB_sites(b_).id        = b_;
                    eNodeB_sites(b_).pos       = [tmp(index) y(index)];
                    eNodeB_sites(b_).parent_centre_pos = micro_centre;
                    eNodeB_sites(b_).site_type = 'micro';
                    eNodeB_sites(b_).sector_centre = sector_centre(s_);
                    eNodeB_sites(b_).O_InH = O_InH;
                    b_ = b_+1;
                    
                    micro_centre = [sector_centre(s_).pos(1)+(ISD/3)/2 sector_centre(s_).pos(2)];
                    tmp = micro_centre(1) + UMi_r * cos(theta);
                    y = micro_centre(2) + UMi_r * sin(theta);
                    index = randperm(361,1);
                    eNodeB_sites(b_)           = network_elements.eNodeB_site;
                    eNodeB_sites(b_).id        = b_;
                    eNodeB_sites(b_).pos       = [tmp(index) y(index)];
                    eNodeB_sites(b_).parent_centre_pos = micro_centre;
                    eNodeB_sites(b_).site_type = 'micro';
                    eNodeB_sites(b_).sector_centre = sector_centre(s_);
                    eNodeB_sites(b_).O_InH = O_InH;
                    b_ = b_+1;
                end
        end
        current_eNodeBs2 = length(eNodeB_sites);
        num_hetnet_sites = current_eNodeBs2 - current_eNodeBs;
end












