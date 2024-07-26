function NR_plot_couplingloss_CDF(SYS_config,UEs,beam_antenna_gain,beam_ue_gain,activate_UE_id)
% ͳ��CL��cdf
% ���������
% SYS_config��LTE����ϵͳ
% UEs��UE�б�
% beam_antenna_gain����վ�˲������ξ���
% beam_ue_gain��UE�˲������ξ���
% activate_UE_id��ͳ�Ƶ�UE��id

for u_=1:length(activate_UE_id)
    macroscopic_pathloss(u_) = UEs(activate_UE_id(u_)).downlink_channel.macroscopic_pathloss;%�Ѿ�����BS_element_gain
    macroscopic_pathloss_900(u_) = UEs(activate_UE_id(u_)).downlink_channel.macroscopic_pathloss_900; % ֻ�����������
    shadow_fading_loss(u_) = UEs(activate_UE_id(u_)).downlink_channel.shadow_fading_pathloss; % ��Ӱ˥��
    penetration_loss(u_) = UEs(activate_UE_id(u_)).downlink_channel.macroscopic_penetration_loss; % ��͸���


    %����BS element gain�Ѿ�����·����
    if SYS_config.antenna_mode~=0
        BS_beamforming_gain(u_)=beam_antenna_gain(activate_UE_id(u_),UEs(activate_UE_id(u_)).attached_eNodeB.eNodeB_id );
    else
        BS_beamforming_gain(u_)=0;
    end

    if SYS_config.antenna_mode==2
        UE_beamforming_gain(u_)=beam_ue_gain(activate_UE_id(u_),UEs(activate_UE_id(u_)).attached_eNodeB.eNodeB_id);
    else
        UE_beamforming_gain(u_)=0;
    end
end
%˵����
%���۵�attach_mode������ã�R1����������+��Ӱ˥��+��͸��ģ�pathloss��
%��attach_mode=1ʱ��R2���pathloss + BS element gain��R3���pathloss + BS array_gain;
%��attach_mode=2ʱ��R2���pathloss + UE element gain, R4���pathloss + UE array_gain
%��attach_mode=3ʱ��R5���pathloss+ UE array gain + BS array gain

pathloss = macroscopic_pathloss_900+shadow_fading_loss+penetration_loss;% R1
pathloss_BS_element_gain = macroscopic_pathloss+shadow_fading_loss+penetration_loss;%BS element gain�Ѿ�����·����(attach_mode=0,1,2,3)
pathloss_BS_antenna_array_gain = pathloss_BS_element_gain - BS_beamforming_gain;%���԰�BS ARRAY GAIN(attach_mode=1)
pathloss_UE_element_gain = pathloss_BS_element_gain - UE_beamforming_gain;%���԰�UE ARRAY GAIN(attach_mode =2)
pathloss_UE_antenna_array_gain = pathloss_BS_element_gain - UE_beamforming_gain-BS_beamforming_gain;%���԰�UE+BS array gain ��attach mode=3��
% pathloss_BS_UE_antenna_array_gain = pathloss_BS_antenna_array_gain - UE_beamforming_gain;

%% Pathloss
pathloss = -pathloss;
min_ = floor(min(pathloss));
max_ = ceil(max(pathloss));
range = min_:1:max_;
h = hist(pathloss,range);
pathloss_CDF = cumsum(h)/sum(h);
if SYS_config.photo
    figure;
    plot(range,pathloss_CDF);
    xlabel('Pathloss��dB��');
    ylabel('CDF');
end

pathloss = sort(pathloss);
index = 1:length(pathloss)/100:length(pathloss);
index = round(index);
R1 = pathloss(index);
R1 = [R1 max(pathloss)];
fid=fopen('.\calibration\R1.txt','wt');
fprintf(fid,'%f\n',R1);
fclose(fid);

%% Pathloss + BS antenna element gain
pathloss_BS_element_gain = -pathloss_BS_element_gain;
min_ = floor(min(pathloss_BS_element_gain));
max_ = ceil(max(pathloss_BS_element_gain));
range = min_:1:max_;
h = hist(pathloss_BS_element_gain,range);
pathloss_element_gain_CDF = cumsum(h)/sum(h);
if SYS_config.photo
    figure;
    plot(range,pathloss_element_gain_CDF);
    xlabel('Pathloss + BS antenna element gain��dB��');
    ylabel('CDF');
end

pathloss_BS_element_gain = sort(pathloss_BS_element_gain);
index = 1:length(pathloss_BS_element_gain)/100:length(pathloss_BS_element_gain);
index = round(index);
R2 = pathloss_BS_element_gain(index);
R2 = [R2 max(pathloss_BS_element_gain)];
fid=fopen('.\calibration\R2.txt','wt');
fprintf(fid,'%f\n',R2);
fclose(fid);

%% Pathloss + BS array gain
pathloss_BS_antenna_array_gain = -pathloss_BS_antenna_array_gain;
min_ = floor(min(pathloss_BS_antenna_array_gain));
max_ = ceil(max(pathloss_BS_antenna_array_gain));
range = min_:1:max_;
h = hist(pathloss_BS_antenna_array_gain,range);
pathloss_array_gain_CDF = cumsum(h)/sum(h);
if SYS_config.photo
    figure;
    plot(range,pathloss_array_gain_CDF);
    xlabel('Pathloss + BS antenna array gain��dB��');
    ylabel('CDF');
end

pathloss_BS_antenna_array_gain = sort(pathloss_BS_antenna_array_gain);
index = 1:length(pathloss_BS_antenna_array_gain)/100:length(pathloss_BS_antenna_array_gain);
index = round(index);
R3 = pathloss_BS_antenna_array_gain(index);
R3 = [R3 max(pathloss_BS_antenna_array_gain)];
fid=fopen('.\calibration\R3.txt','wt');
fprintf(fid,'%f\n',R3);
fclose(fid);

%% Pathloss + UE antenna element gain
pathloss_UE_element_gain = -pathloss_UE_element_gain;
min_ = floor(min(pathloss_UE_element_gain));
max_ = ceil(max(pathloss_UE_element_gain));
range = min_:1:max_;
h = hist(pathloss_UE_element_gain,range);
tmp_CDF = cumsum(h)/sum(h);
if SYS_config.photo
    figure;
    plot(range,tmp_CDF);
    xlabel('Pathloss + UE antenna element gain��dB��');
    ylabel('CDF');
end

pathloss_UE_element_gain = sort(pathloss_UE_element_gain);
index = 1:length(pathloss_UE_element_gain)/100:length(pathloss_UE_element_gain);
index = round(index);
R4 = pathloss_UE_element_gain(index);
R4 = [R4 max(pathloss_UE_element_gain)];
fid=fopen('.\calibration\R4.txt','wt');
fprintf(fid,'%f\n',R4);
fclose(fid);

%% Pathloss + UE array gain
pathloss_UE_antenna_array_gain = -pathloss_UE_antenna_array_gain;
min_ = floor(min(pathloss_UE_antenna_array_gain));
max_ = ceil(max(pathloss_UE_antenna_array_gain));
range = min_:1:max_;
h = hist(pathloss_UE_antenna_array_gain,range);
tmp_CDF = cumsum(h)/sum(h);
if SYS_config.photo
    figure;
    plot(range,tmp_CDF);
    xlabel('Pathloss + UE antenna array gain��dB��');
    ylabel('CDF');
end

pathloss_UE_antenna_array_gain = sort(pathloss_UE_antenna_array_gain);
index = 1:length(pathloss_UE_antenna_array_gain)/100:length(pathloss_UE_antenna_array_gain);
index = round(index);
R5 = pathloss_UE_antenna_array_gain(index);
R5 = [R5 max(pathloss_UE_antenna_array_gain)];
fid=fopen('.\calibration\R5.txt','wt');
fprintf(fid,'%f\n',R5);
fclose(fid);






