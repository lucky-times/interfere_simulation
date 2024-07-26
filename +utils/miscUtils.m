classdef miscUtils
    % ����ת���Ƕȵ�һЩС����
    % �����ƽ̨ԭ����һ�������������Ҳ�Ž�����
    
    properties
    end
    
    methods(Static)
        
        % ���Ƕ�ת��Ϊ[0 359]
        function mod_angle = wrapTo359(the_angle)
            mod_angle = mod(the_angle, 360);
        end
        
        % ���Ƕ�ת��Ϊ[0 180]
        function mod_angle = wrapTo180(the_angle)
            the_angle(the_angle-180>0) = the_angle(the_angle-180>0)-360;
            mod_angle = the_angle;
        end
        
        % ���Ƕ�ת��Ϊ[-180 180]
        function angle = wrapToAll180(angle_ini)
            % ���ڽ��Ƕ�ת��Ϊ-pi~pi
            % ���������
            % angle_ini��ת��ǰ�ĽǶȣ��ȣ�
            %
            % ���������
            % angle��ת����ĽǶȣ��ȣ�
            
            angle_ini=angle_ini+180;
            angle_ini=mod(angle_ini,360);
            angle=angle_ini-180;
        end
        
        % ����element���е�clear����
        function tidy_up_memory_before_closing(UEs,eNodeBs_sectors,eNodeBs)
            for u_=1:length(UEs)
                if ~isstruct(UEs(u_))
                    UEs(u_).clear;
                end
            end
            
            for bs_=1:length(eNodeBs_sectors)
                eNodeBs_sectors(bs_).clear;
            end
            
            for site_=1:length(eNodeBs)
                eNodeBs(site_).clear;
            end
        end
        
        % ��GUI�и���listbox��ѡ���С��ѡ��������ЩС����UE
        function [UE_ids cell_sum_throughput] = get_UEs_in_given_cells(cell_ids,the_UE_traces)
            N_UEs  = length(the_UE_traces);
            cell_sum_throughput = zeros(1,length(cell_ids));
            UE_ids = false(1,N_UEs);
            for c_idx=1:length(cell_ids)
                c_ = cell_ids(c_idx);
                for u_=1:N_UEs
                    if ~isempty(find(the_UE_traces(u_).attached_eNodeB==c_,1,'first'))
                        UE_ids(u_) = true;
                    end
                end
            end
        end

        %ȫ������ϵת��Ϊ��������ϵ
        function [phi2,theta2]=global2local(phi1,theta1,alpha,beta,gamma)%alpha��ʾ��z����ת�Ƕȣ�beta��ʾy����ת�Ƕȣ�gamma��ʾx����ת�Ƕ�,phi��ʾˮƽ�ǣ�theta��ʾ�����
            
            R1=[cosd(alpha),sind(alpha),0;
                -sind(alpha),cosd(alpha),0;
                0,0,1];%��z����ת����
            R2=[cosd(beta),0,-sind(beta);
                0,1,0;
                sind(beta),0,cosd(beta)];%��y����ת����
            R3=[1,0,0;
                0,cosd(gamma),sind(gamma);
                0,-sind(gamma),cosd(gamma)];%��x����ת����
            R = R1*R2*R3;
            direc=[sind(theta1).*cosd(phi1);
                sind(theta1).*sind(phi1);
                cosd(theta1)];
            tran=inv(R)*direc;
            theta2=acosd([0,0,1]*tran);%���ؽǶ�ֵ
            phi2=angle([1,1j,0]*tran)/pi*180;
        end

    end
end

