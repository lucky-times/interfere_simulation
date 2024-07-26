function NR_throughput_asy(SYS_config,ul_dl_wideband_loop_SINR,dl_ul_wideband_loop_SINR,activate_UE_id)
% ����ӳ��SINR��������(�����в�ͬ��ͳ��ģ��)
% ���������
% SYS_config��LTE����ϵͳ
% UEs��UE�б�
% wideband_loop_SINR�����и�������SINR

n_UE_served_per_BS = SYS_config.n_UE_served_per_BS;
n_UE = length(activate_UE_id);
bandwidth = SYS_config.bandwidth/n_UE_served_per_BS;
n_snapshot = SYS_config.UE_per_eNodeB/n_UE_served_per_BS;
[loop_size,~,~] = size(ul_dl_wideband_loop_SINR);

%����36942��ӳ�乫ʽʵ��ӳ��
%% ���и�������
DL_SINR = ul_dl_wideband_loop_SINR;
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

%% ���и�������
UL_SINR = dl_ul_wideband_loop_SINR;
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
fid=fopen('.\calibration\Throughput asynchronization.txt','wt');
for i=1:loop_size
    if i==loop_size
        str = 'no interference';
    else
        str = num2str(SYS_config.ACIR_lower+SYS_config.loop_step*(i-1));
    end
    msg = str;
    fprintf(fid,'%s\n',msg);
    msg = ['UE to UE average per cell: ',num2str(ave_cell_DL_thput(i,:)),' ','UE to UE edge per cell: ',num2str(cell_DL_thput_5(i,:))];
    fprintf(fid,'%s\n',msg);
    msg = ['BS to BS average per cell: ',num2str(ave_cell_UL_thput(i,:)),' ','BS to BS edge per cell: ',num2str(cell_UL_thput_5(i,:))];
    fprintf(fid,'%s\n',msg);
end
fclose(fid);

fid=fopen('.\calibration\Throughput loss asynchronization.txt','wt');
percentage = zeros(4,loop_size-1);

for i=1:(loop_size-1)
    percentage(1,i) =  (ave_cell_DL_thput(loop_size,:) - ave_cell_DL_thput(i,:))/ave_cell_DL_thput(loop_size,:) * 100;
    percentage(2,i) =  (cell_DL_thput_5(loop_size,:) - cell_DL_thput_5(i,:))/cell_DL_thput_5(loop_size,:) * 100;
    percentage(3,i) =  (ave_cell_UL_thput(loop_size,:) - ave_cell_UL_thput(i,:))/ave_cell_UL_thput(loop_size,:) * 100;
    percentage(4,i) =  (cell_UL_thput_5(loop_size,:) - cell_UL_thput_5(i,:))/cell_UL_thput_5(loop_size,:) * 100;
end
fprintf(fid,'%s\n','UE to UE average per cell');
fprintf(fid,'%f ',percentage(1,:));
fprintf(fid,'\n');
fprintf(fid,'%s\n','UE to UE edge per cell');
fprintf(fid,'%f ',percentage(2,:));
fprintf(fid,'\n');
fprintf(fid,'%s\n','BS to BS average per cell');
fprintf(fid,'%f ',percentage(3,:));
fprintf(fid,'\n');
fprintf(fid,'%s\n','BS to BS edge per cell');
fprintf(fid,'%f ',percentage(4,:));
fprintf(fid,'\n');
fclose(fid);




ACIR=(SYS_config.ACIR_lower+SYS_config.loop_step.*[0:1:loop_size-2]);%-30;
%X=(SYS_config.X_lower-30):SYS_config.X_loop_step:(SYS_config.X_upper-30);
%X=(SYS_config.X_lower):SYS_config.X_loop_step:(SYS_config.X_upper);
%X=ACIR-30;
figure;
plot(ACIR,percentage(1,:),'-o','color','r');
hold on
plot(ACIR,percentage(2,:),'-s','color','b');
%axis([-inf,inf,0,1]);
%hold on
%y=5;
%plot(ACIR,y)
grid on;
xlabel('ACIR��dB��');
ylabel('Throughput Loss (%)');
title('UL to DL Throughput Loss');
legend('UE to UE average throughput loss','UE to UE edge throughput loss')
figure;
plot(ACIR,percentage(3,:),'-o','color','r');
hold on
plot(ACIR,percentage(4,:),'-s','color','b');
%axis([-inf,inf,0,1]);
grid on;
xlabel('ACIR��dB��');
ylabel('Throughput Loss (%)');
title('DL to UL Throughput Loss');
legend('BS to BS average throughput loss','BS to BS cell edge throughput loss')



