function NR_throughput(SYS_config,wideband_loop_SINR,ul_wideband_loop_SINR,activate_UE_id)
% ����ӳ��SINR��������
% ���������
% SYS_config��LTE����ϵͳ
% UEs��UE�б�
% wideband_loop_SINR������SINR
% ul_wideband_loop_SINR�� ����SINR

n_UE_served_per_BS = SYS_config.n_UE_served_per_BS;
n_UE = length(activate_UE_id);
bandwidth = SYS_config.bandwidth/n_UE_served_per_BS;
n_snapshot = SYS_config.UE_per_eNodeB/n_UE_served_per_BS;
[loop_size,~,~] = size(wideband_loop_SINR);

%����36942��ӳ�乫ʽʵ��ӳ��
%% ����
DL_SINR = wideband_loop_SINR;
for u_=1:n_UE
    for i=1:loop_size
        if DL_SINR(i,:,u_)<-10 % ��SINRС��-10dB��������Ϊ0
            DL_thput(i,:,u_) = 0;
        elseif DL_SINR(i,:,u_)>30 % ��SINR����30dB��������Ϊ���������
            DL_thput(i,:,u_) = 0.6 * bandwidth/1e6*log2(1+10^(0.1*30));
        else % ��SINR��-10~30֮�䣬ͨ����ũ��ʽ�õ�������
            DL_thput(i,:,u_) = 0.6 * bandwidth/1e6*log2(1+10^(0.1*DL_SINR(i,:,u_)));
        end
    end
end
% ƽ��������
ave_DL_thput = sum(DL_thput,3)/n_UE;
% ��Ե������
if n_snapshot~=1
    % ���ӣ���һ�ο��ջ�վ1����UE1��2����վ2����UE11��12������Ȼ��ͳ�Ƹÿ����±�ԵUE������������������п����±�ԵUE������������ƽ��
    for j=1:n_snapshot
        index_UE = [];
        for a_ = 1:n_UE_served_per_BS
            if isempty(index_UE)
                index_UE = (j*n_UE_served_per_BS)-(n_UE_served_per_BS-a_):SYS_config.UE_per_eNodeB:n_UE;
            else
                index_UE = [index_UE,(j*n_UE_served_per_BS)-(n_UE_served_per_BS-a_):SYS_config.UE_per_eNodeB:n_UE ];
            end
        end
        DL_thput_divi(:,:,:,j) = DL_thput(:,:,index_UE);
    end
    DL_thput_divi_sort = sort(DL_thput_divi,3);
    
    index_5 = ceil(n_UE/n_snapshot*0.05);
    DL_thput_divi_ave = sum(DL_thput_divi_sort,4)/n_snapshot;
    
    DL_thput_5 = DL_thput_divi_ave(:,:,index_5);
else
    DL_thput_sort = sort(DL_thput,3);
    index_5 = ceil(n_UE*0.05);
    DL_thput_5 = DL_thput_sort(:,:,index_5);
end

%% ����
UL_SINR = ul_wideband_loop_SINR;
for u_=1:n_UE
    for i=1:loop_size
        if UL_SINR(i,:,u_)<-10
            UL_thput(i,:,u_) = 0;
        elseif UL_SINR(i,:,u_)>22
            UL_thput(i,:,u_) = 0.4 * bandwidth/1e6*log2(1+10^(0.1*22));
        else
            UL_thput(i,:,u_) = 0.4 * bandwidth/1e6*log2(1+10^(0.1*UL_SINR(i,:,u_)));
        end
    end
end

ave_UL_thput = sum(UL_thput,3)/n_UE;

if n_snapshot~=1
    for j=1:n_snapshot
        index_UE = [];
        for a_ = 1:n_UE_served_per_BS
            if isempty(index_UE)
                index_UE = (j*n_UE_served_per_BS)-(n_UE_served_per_BS-a_):SYS_config.UE_per_eNodeB:n_UE;
            else
                index_UE = [index_UE,(j*n_UE_served_per_BS)-(n_UE_served_per_BS-a_):SYS_config.UE_per_eNodeB:n_UE ];
            end
        end
        UL_thput_divi(:,:,:,j) = UL_thput(:,:,index_UE);
    end
    UL_thput_divi_sort = sort(UL_thput_divi,3);
    
    index_5 = ceil(n_UE/n_snapshot*0.05);
    UL_thput_divi_ave = sum(UL_thput_divi_sort,4)/n_snapshot;
    
    UL_thput_5 = UL_thput_divi_ave(:,:,index_5);
else
    UL_thput_sort = sort(UL_thput,3);
    index_5 = ceil(n_UE*0.05);
    UL_thput_5 = UL_thput_sort(:,:,index_5);
end

%% ͳ��
ave_cell_DL_thput = ave_DL_thput*n_UE_served_per_BS;
cell_DL_thput_5 = DL_thput_5*n_UE_served_per_BS;
ave_cell_UL_thput = ave_UL_thput*n_UE_served_per_BS;
cell_UL_thput_5 = UL_thput_5*n_UE_served_per_BS;
%% ��������ʧ
fid=fopen('.\calibration\Throughput.txt','wt');
for i=1:loop_size
    if i==loop_size
        str = 'no interference';
    else
        str = num2str(5+SYS_config.loop_step*(i-1));
    end
    msg = str;
    fprintf(fid,'%s\n',msg);
    msg = ['DL average per cell: ',num2str(ave_cell_DL_thput(i,:)),' ','DL edge per cell: ',num2str(cell_DL_thput_5(i,:))];
    fprintf(fid,'%s\n',msg);
    msg = ['UL average per cell: ',num2str(ave_cell_UL_thput(i,:)),' ','UL edge per cell: ',num2str(cell_UL_thput_5(i,:))];
    fprintf(fid,'%s\n',msg);
end
fclose(fid);

fid=fopen('.\calibration\Throughput loss.txt','wt');
percentage = zeros(4,loop_size-1);

for i=1:(loop_size-1)
    percentage(1,i) =  (ave_cell_DL_thput(loop_size,:) - ave_cell_DL_thput(i,:))/ave_cell_DL_thput(loop_size,:) * 100;
    percentage(2,i) =  (cell_DL_thput_5(loop_size,:) - cell_DL_thput_5(i,:))/cell_DL_thput_5(loop_size,:) * 100;
    percentage(3,i) =  (ave_cell_UL_thput(loop_size,:) - ave_cell_UL_thput(i,:))/ave_cell_UL_thput(loop_size,:) * 100;
    percentage(4,i) =  (cell_UL_thput_5(loop_size,:) - cell_UL_thput_5(i,:))/cell_UL_thput_5(loop_size,:) * 100;
end
fprintf(fid,'%s\n','DL average');
fprintf(fid,'%f ',percentage(1,:));
fprintf(fid,'\n');
fprintf(fid,'%s\n','DL edge');
fprintf(fid,'%f ',percentage(2,:));
fprintf(fid,'\n');
fprintf(fid,'%s\n','UL average');
fprintf(fid,'%f ',percentage(3,:));
fprintf(fid,'\n');
fprintf(fid,'%s\n','UL edge');
fprintf(fid,'%f ',percentage(4,:));
fprintf(fid,'\n');
fclose(fid);

%% �����ʧ
user_loss = zeros(2,(loop_size - 1));
% ����
m = find(DL_thput(loop_size,:,:) == 0);
user_connect_noacir = 1 - length(m)/n_UE;
% ����
m = find(UL_thput(loop_size,:,:) == 0);
user_connect_noacir_ul = 1 - length(m)/n_UE;
for i=1:(loop_size - 1)
   m = find(DL_thput(i,:,:) == 0);
   user_loss(1,i) = (user_connect_noacir - (1 - length(m)/n_UE))*100;
   m = find(UL_thput(i,:,:) == 0);
   user_loss(2,i) = (user_connect_noacir_ul - (1 - length(m)/n_UE))*100;
end
    
fid=fopen('.\calibration\User loss.txt','wt');
fprintf(fid,'%s\n','DL');
fprintf(fid,'%f ',user_loss(1,:));
fprintf(fid,'\n');
fprintf(fid,'%s\n','UL');
fprintf(fid,'%f ',user_loss(2,:));
fclose(fid);








