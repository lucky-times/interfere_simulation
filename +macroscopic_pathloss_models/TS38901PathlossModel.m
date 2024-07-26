classdef TS38901PathlossModel < macroscopic_pathloss_models.generalPathlossModel
    % ͨ��38900����·����ģ�������ģ�
    properties
        frequency     % ϵͳƵ��
        environment   % ����
        h_ut          %ut�����߸߶�
        h_bs          %base�����߸߶�
        pathloss_function_handle

        PL_UMa_LOS    %UMA LOS�ֲ����
        PL_UMi_LOS    %UMI LOS�ֲ����
        PL_InH_LOS    %INH LOS�ֲ����
        PL_RMa_LOS    %RMA LOS�ֲ����
    end

    methods
        function obj = TS38901PathlossModel(frequency,environment,SYS_config) %��ʼ��
            obj.frequency = frequency;
            obj.environment = environment;
            obj.name = 'TS 38.901';
            switch environment %ѡ�񳡾���������Ӧ�Ĳ���
                case 'urban_macro'
                    obj.name = [obj.name ' urban macro'];
                    obj.h_ut=SYS_config.UE_height;
                    obj.h_bs=25;
                    obj.pathloss_function_handle = 1;
                case 'urban_micro'
                    obj.name = [obj.name ' urban micro'];
                    obj.h_ut=1.5;
                    obj.h_bs=10;
                    obj.pathloss_function_handle = 2;
                case 'indoor'
                    obj.name = [obj.name ' indoor'];
                    obj.h_ut=1;% �̶�
                    obj.h_bs=3;% �̶�
                    obj.pathloss_function_handle = 3;
                case 'rural_macro'
                    obj.name = [obj.name ' rural macro'];
                    obj.h_ut=1.5;
                    obj.h_bs=35;
                    obj.pathloss_function_handle = 4;
            end
        end

        function [pathloss_in_db,comp] = pathloss(obj,d_2d,h_ut,d_2d_in)
            %ȷ������������LOS�ֲ��������㴫�����
            switch obj.pathloss_function_handle
                case 1
                    [pathloss_in_db,comp] = obj.pathloss_urban_macro_LOS_NLOS(d_2d,h_ut,d_2d_in);
                case 2
                    [pathloss_in_db,comp] = obj.pathloss_urban_micro_LOS_NLOS(d_2d,h_ut,d_2d_in);
                case 3
                    [pathloss_in_db,comp] = obj.pathloss_indoor_LOS_NLOS(d_2d,h_ut);
                case 4
                    [pathloss_in_db,comp] = obj.pathloss_rural_macro_LOS_NLOS(d_2d,h_ut,d_2d_in);
            end
        end
        %% ����RMA����NLOS��ʽ������Ӧ�������
        % ����ֱ�Ϊ��ά������󣬸߶Ⱦ����Լ�������Ծ������
        function pl = pathloss_rural_macro_NLOS(obj,d_2d,h_ut,d_2d_in)
            d_3d=sqrt(d_2d.^2+(obj.h_bs-h_ut).^2);
            d_3d_in = d_3d.*d_2d_in./d_2d;
            d_3d_out = d_3d-d_3d_in;

            w = 20;
            h = 5;
            pl_NLOS=161.04-7.1*log10(w)+7.5*log10(h)-(24.37-3.7*((h/obj.h_bs)^2))*log10(obj.h_bs)+(43.42-3.1*log10(obj.h_bs)).*(log10(d_3d_out)-3)+20*log10(obj.frequency/1e9)-(3.2.*(log10(11.75.*h_ut).^2-4.97));
            pl_LOS = obj.PL_RMa_LOS;
            pl=max(pl_LOS,pl_NLOS);
        end

        function [pl,comp2] = pathloss_rural_macro(obj,d_2d,h_ut,d_2d_in)
            d_bp = 2*pi*obj.h_bs.*h_ut*obj.frequency/3e8;
            d_2d_out = d_2d - d_2d_in;
            d_3d=sqrt(d_2d.^2+(obj.h_bs-h_ut).^2);
            d_3d_in = d_3d.*d_2d_in./d_2d;
            d_3d_out = d_3d - d_3d_in;

            pl = zeros(size(d_2d));
            comp2 = zeros(size(d_2d));
            [row,col] = size(d_2d);
            h = 5;
            PL1 = @(d) 20*log10(40*pi*d*(obj.frequency/1e9)/3)+min(0.03*(h^1.72),10)*log10(d)-min(0.044*(h^1.72),14.77)+0.002*log10(h)*d;
            for i = 1:row*col
                if d_2d_out(i) <= d_bp(i) && d_2d_out(i) >= 10
                    pl(i) = PL1(d_3d_out(i));
                    comp2(i) = -1;
                elseif d_2d_out(i) <= 10000 && d_2d_out(i) > d_bp(i)
                    pl(i) = PL1(d_bp(i)) + 40*log10(d_3d_out(i)/d_bp(i));
                end
            end
            obj.PL_RMa_LOS = pl;
        end

        function [pl,comp]=pathloss_rural_macro_LOS_NLOS(obj,d_2d,h_ut,d_2d_in)
            [pl_LOS,comp2]=obj.pathloss_rural_macro(d_2d,h_ut,d_2d_in);
            pl_NLOS=obj.pathloss_rural_macro_NLOS(d_2d,h_ut,d_2d_in);
            [row,col] = size(d_2d);
            d_2d_out = d_2d - d_2d_in;

            LOS_posibility = ones(size(d_2d));
            [m,n] = find(d_2d_out>10);
            for i = 1:length(m)
                LOS_posibility(m(i),n(i)) =  exp(-(d_2d_out(m(i),n(i))-10)/1000);% LOS����
            end

            r_rnd=rand(row,col);
            compare=r_rnd>LOS_posibility;%���compare=1�����NLOS�����compare=0������LOS;
            [NLOS_x,NLOS_y]=find(compare==1);
            NLOS_length=length(NLOS_x);
            for ii=1:NLOS_length
                pl(NLOS_x(ii),NLOS_y(ii))=pl_NLOS(NLOS_x(ii),NLOS_y(ii));
            end
            [LOS_x,LOS_y]=find(compare==0);
            LOS_length=length(LOS_x);
            for ij=1:LOS_length
                pl(LOS_x(ij),LOS_y(ij))=pl_LOS(LOS_x(ij),LOS_y(ij));
            end
            [m,n] = find(compare==0);
            for i = 1:length(m)
                compare(m(i),n(i)) = comp2(m(i),n(i));
            end
            comp=compare;
        end

        %% ����UMA����NLOS��ʽ������Ӧ�������
        % ����ֱ�Ϊ��ά������󣬸߶Ⱦ����Լ�������Ծ������
        function pl =pathloss_urban_macro_NLOS(obj,d_2d,h_ut,d_2d_in)
            d_3d=sqrt(d_2d.^2+(obj.h_bs-h_ut).^2);
            d_3d_in = d_3d.*d_2d_in./d_2d;
            d_3d_out = d_3d-d_3d_in;
            pl_NLOS=13.54+39.08*log10(d_3d_out)+20*log10(obj.frequency/1e9)-0.6*(h_ut-1.5);
            pl_LOS = obj.PL_UMa_LOS;
            pl=max(pl_LOS,pl_NLOS);
        end

        function pl = pathloss_urban_macro(obj,d_2d,h_ut,d_2d_in)

            [row,col] = size(d_2d);
            d_2d_out = d_2d - d_2d_in; %�������ص�ര�ľ���
            g = zeros(size(h_ut));
            C = zeros(size(h_ut));
            for i = 1:row*col
                if d_2d_in(i) == 0
                    if d_2d(i)<=18
                        g(i) = 0;
                    else
                        g(i) = 5/4*(d_2d(i)/100)^3*exp(-d_2d(i)/150);
                    end
                else
                    if d_2d_out(i)<=18
                        g(i) = 0;
                    else
                        g(i) = 5/4*(d_2d(i)/100)^3*exp(-d_2d(i)/150);
                    end
                end
            end
            [m,n] = find(h_ut(:,:)>=13);
            for i=1:length(m)
                C(m(i),n(i)) = (((h_ut(m(i),n(i))-13)/10)^1.5)*g(m(i),n(i));
            end
            p = 1./(1+C);
            h_e = zeros(size(h_ut));
            pl = zeros(size(h_ut));
            tmp = rand(row,col);
            for i=1:row*col
                if h_ut(i)>=13
                    opt = 12:3:h_ut(i)-1;
                    select = randperm(length(opt));
                    h_e(i) = opt(select(1));
                else
                    h_e(i) = 1;
                end
            end
            h_e(find(tmp<=p)) = 1;
            h_bs_1 = obj.h_bs - h_e;
            h_ut_1 = h_ut-h_e;
            d_bp_1 = 4*h_bs_1.*h_ut_1*obj.frequency/3e8; %����ֽ��

            d_3d=sqrt(d_2d.^2+(obj.h_bs-h_ut).^2);
            d_3d_in = d_3d.*d_2d_in./d_2d;
            d_3d_out = d_3d - d_3d_in;
            for i=1:row*col
                if d_bp_1(i)<10
                    error('it is too near');
                else
                    if d_2d_out(i)<d_bp_1(i)
                        pl(i) = 28 + 22*log10(d_3d_out(i)) + 20*log10( obj.frequency/1e9);
                    else
                        pl(i) = 28 + 40*log10(d_3d_out(i)) + 20*log10( obj.frequency/1e9)-9*log10(d_bp_1(i)^2+(obj.h_bs-h_ut(i))^2);
                    end
                end
            end
            obj.PL_UMa_LOS = pl;
        end
        %% ����UMAȷ��LOS��NLOS�ı������������ϵõ�����·��
        % ��������ֱ�Ϊ�����ص㵽��С����·��ֲ�����ͼ��LOS�ֲ�����
        function [pl,comp]=pathloss_urban_macro_LOS_NLOS(obj,d_2d,h_ut,d_2d_in)
            %����UAV������ȫ������los��������nlos���
            pl_LOS=obj.pathloss_urban_macro(d_2d,h_ut,d_2d_in);
            pl_NLOS=obj.pathloss_urban_macro_NLOS(d_2d,h_ut,d_2d_in);
            [row,col] = size(d_2d);
            d_2d_out = d_2d - d_2d_in;
            C = zeros(size(h_ut));
            [m,n] = find(h_ut(:,:)>=13);
            for i=1:length(m)
                C(m(i),n(i)) = ((h_ut(m(i),n(i))-13)/10)^1.5;
            end

            % ����LOS�ĸ���
            LOS_posibility = ones(size(d_2d_out));
            [m,n] = find(d_2d_out>18);
            for i = 1:length(m)
                LOS_posibility(m(i),n(i)) = (18/d_2d_out(m(i),n(i))+exp(-d_2d_out(m(i),n(i))/63)*(1-18/d_2d_out(m(i),n(i))))*(1+C(m(i),n(i))*5/4*(d_2d_out(m(i),n(i))/100)^3*exp(-d_2d_out(m(i),n(i))/150));
            end

            r_rnd=rand(row,col);
            %             compare=r_rnd>LOS_posibility;%���compare=1�����NLOS�����compare=0������LOS;
            compare = zeros(row, col);
            [NLOS_x,NLOS_y]=find(compare==1);
            NLOS_length=length(NLOS_x);
            for ii=1:NLOS_length
                pl(NLOS_x(ii),NLOS_y(ii))=pl_NLOS(NLOS_x(ii),NLOS_y(ii));
            end
            [LOS_x,LOS_y]=find(compare==0);
            LOS_length=length(LOS_x);
            for ij=1:LOS_length
                pl(LOS_x(ij),LOS_y(ij))=pl_LOS(LOS_x(ij),LOS_y(ij));
            end
            comp=compare;
        end
        %% ����UMI����NLOS��ʽ������Ӧ�������
        function pl =pathloss_urban_micro_NLOS(obj,d_2d,h_ut,d_2d_in)
            d_3d=sqrt(d_2d.^2+(obj.h_bs-h_ut).^2);
            d_3d_in = d_3d.*d_2d_in./d_2d;
            d_3d_out = d_3d-d_3d_in;
            pl_NLOS=35.3*log10(d_3d_out)+22.4+21.3*log10(obj.frequency/1e9)-0.3*(h_ut-1.5);
            pl_LOS = obj.PL_UMi_LOS;
            pl=max(pl_LOS,pl_NLOS);
        end

        function pl = pathloss_urban_micro(obj,d_2d,h_ut,d_2d_in)
            [row,col] = size(h_ut);
            pl = zeros(size(h_ut));

            h_bs_1 = obj.h_bs -1;
            h_ut_1 = h_ut-1;
            d_bp_1 = 4*h_bs_1.*h_ut_1*obj.frequency/3e8;

            d_3d=sqrt(d_2d.^2+(obj.h_bs-obj.h_ut)^2);
            d_3d_in = d_3d.*d_2d_in./d_2d;
            d_3d_out = d_3d-d_3d_in;
            d_2d_out = d_2d-d_2d_in;

            for i=1:row*col
                if d_bp_1(i)<10
                    error('it is too near');
                else
                    if d_2d_out(i)<d_bp_1(i)
                        pl(i) = 32.4 + 21*log10(d_3d_out(i)) + 20*log10( obj.frequency/1e9);
                    else
                        pl(i) = 32.4 + 40*log10(d_3d_out(i)) + 20*log10( obj.frequency/1e9)-9.5*log10(d_bp_1(i)^2+(obj.h_bs-h_ut(i))^2);
                    end
                end
            end
            obj.PL_UMi_LOS = pl;
        end
        %% ����UMiȷ��LOS��NLOS�ı������������ϵõ�����·��
        function [pl,comp]=pathloss_urban_micro_LOS_NLOS(obj,d_2d,h_ut,d_2d_in)
            [row,col] = size(h_ut);
            pl_LOS=obj.pathloss_urban_micro(d_2d,h_ut,d_2d_in);
            pl_NLOS=obj.pathloss_urban_micro_NLOS(d_2d,h_ut,d_2d_in);
            d_2d_out = d_2d - d_2d_in;

            % ����LOS�ĸ���
            LOS_posibility = ones(size(d_2d_out));
            [m,n] = find(d_2d_out>18);
            for i = 1:length(m)
                LOS_posibility(m(i),n(i)) = 18/d_2d_out(m(i),n(i))+exp(-d_2d_out(m(i),n(i))/36)*(1-18/d_2d_out(m(i),n(i)));
            end

            r_rnd=rand(row,col);
            compare=r_rnd>LOS_posibility;%���compare=1�����NLOS�����compare=0������LOS;
            [NLOS_x,NLOS_y]=find(compare==1);
            NLOS_length=length(NLOS_x);
            for ii=1:NLOS_length
                pl(NLOS_x(ii),NLOS_y(ii))=pl_NLOS(NLOS_x(ii),NLOS_y(ii));
            end
            [LOS_x,LOS_y]=find(compare==0);
            LOS_length=length(LOS_x);
            for ij=1:LOS_length
                pl(LOS_x(ij),LOS_y(ij))=pl_LOS(LOS_x(ij),LOS_y(ij));
            end
            comp=compare;
        end
        %% ����InH����NLOS��ʽ������Ӧ�������
        function pl = pathloss_indoor_NLOS(obj,d_2d,h_ut)
            d_3d = sqrt((obj.h_bs-h_ut).^2 + d_2d.^2);
            pl_NLOS = 38.3*log10(d_3d) + 17.3 + 24.9*log10(obj.frequency/1e9);
            pl_LOS = obj.PL_InH_LOS;
            pl=max(pl_LOS,pl_NLOS);
        end
        %% ����InH����LOS��ʽ������Ӧ�������
        function pl = pathloss_indoor(obj,d_2d,h_ut)
            d_3d = sqrt((obj.h_bs-h_ut).^2 + d_2d.^2);
            pl = 32.4 + 17.3*log10(d_3d) + 20*log10(obj.frequency/1e9);
            obj.PL_InH_LOS = pl;
        end
        %% ����InHȷ��LOS��NLOS�ı������������ϵõ�����·��
        function [pl,comp]=pathloss_indoor_LOS_NLOS(obj,d_2d,h_ut)
            [row,col] = size(d_2d);
            pl_LOS=obj.pathloss_indoor(d_2d,h_ut);
            pl_NLOS=obj.pathloss_indoor_NLOS(d_2d,h_ut);

            LOS_posibility = zeros(size(d_2d));
            for i = 1:row*col
                if d_2d(i)<=5
                    LOS_posibility(i) = 1;
                elseif 5<d_2d(i) && d_2d(i)<=49
                    LOS_posibility(i) = exp(-(d_2d(i)-5)/70.8);
                else
                    LOS_posibility(i) = exp(-(d_2d(i)-49)/211.7)*0.54;
                end
            end
            r_rnd=rand(row,col);
            compare=r_rnd>LOS_posibility;%���compare=1�����NLOS�����compare=0������LOS;
            [NLOS_x,NLOS_y]=find(compare==1);
            NLOS_length=length(NLOS_x);
            for ii=1:NLOS_length
                pl(NLOS_x(ii),NLOS_y(ii))=pl_NLOS(NLOS_x(ii),NLOS_y(ii));
            end
            [LOS_x,LOS_y]=find(compare==0);
            LOS_length=length(LOS_x);
            for ij=1:LOS_length
                pl(LOS_x(ij),LOS_y(ij))=pl_LOS(LOS_x(ij),LOS_y(ij));
            end
            comp=compare;
        end
    end
end