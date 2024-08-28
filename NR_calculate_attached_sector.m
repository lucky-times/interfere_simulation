function [eNodeB_sector_id] = NR_calculate_attached_sector(SYS_config, eNodeB_sites,eNodeB_sectors, UE)
%UNTITLED3 此处提供此函数的摘要
%   此处提供详细说明
UE_pos = UE.pos;
UE_height = UE.height;

% 根据天线增益判断归属小区
m = length(eNodeB_sites);
n = length(eNodeB_sectors);
BS_pos = zeros(m, 2);
for i = 1:m
    BS_pos(i, :) = eNodeB_sites(i).pos;
end
antenna = antennas.BSAntenna(6, 8, 1, 0);
sector_azimuth = [0;120;240];
d = zeros(m,1);
for i = 1:m
    d(i) = pdist2(eNodeB_sites(i).pos, UE_pos);
end
% 计算垂直天顶角度
D = repmat(d, 1, 3);
h = zeros(m,1);
for i = 1:m
    h(i) = UE_height - eNodeB_sectors((i-1)*3+1).tx_height;
end
H = repmat(h, 1, 3);
el = atand(D./H);
el = el';
% 计算水平方位角度
ang = (180/pi)*(atan2((UE_pos(2)-BS_pos(:,2)),(UE_pos(1)-BS_pos(:,1))));
Ang = repmat(ang', 3, 1)';
sector_azimuth = zeros(size(Ang));
for ii = 1:n
    i = ceil(ii/3);
    j = ii-(i-1)*3;
    sector_azimuth(i,j) = eNodeB_sectors(ii).azimuth;
end
az = Ang - sector_azimuth;%mod(Ang, 360)
horizontal_angle_grid = utils.miscUtils.wrapTo359(az);
horizontal_angle_grid_s = utils.miscUtils.wrapTo180(horizontal_angle_grid);
az = horizontal_angle_grid_s';
% save('abc.mat','az', 'el', 'elementGain')

elementGain = zeros(size(az));
for i = 1:size(az, 1)
    for j = 1:size(az, 2)
        elementGain(i,j) = antenna.elementPattern(el(i,j), az(i,j), 0);
    end
end
% elementGain = antenna.elementPattern(el, az);
[gain,  index] = sort(reshape(elementGain, [numel(elementGain), 1]));
eNodeB_sector_id = index(end);


% % 根据垂直角度和水平角度判断归属小区
% m = length(eNodeB_sites);
% BS_pos = zeros(m, 2);
% for i = 1 : m
%     BS_pos(i, :) = eNodeB_sites(i).pos;
% end
% 
% d = zeros(m,1);
% for i = 1:m
%     d(i) = pdist2(BS_pos(i,:), UE_pos);
% end
% [~,  index] = sort(d);
% 
% eNodeB_order = 0;
% for i = 1:m
%     alpha = atand((UE_height - eNodeB_sectors((i-1)*3+1).tx_height) / d(index(i)));
%     if(alpha <= (SYS_config.tilt+18) && alpha > SYS_config.tilt)
%         eNodeB_order = index(i);
%         break;
%     end
% end
% if eNodeB_order ==0
%     error("归属基站判别出错！")
% end
% 
% attached_sectors = [];
% for i = 1:length(eNodeB_sectors)
%     if eNodeB_sectors(i).parent_eNodeB.id == eNodeB_sites(eNodeB_order).id
%         attached_sectors = [attached_sectors;eNodeB_sectors(i)];
%     end
% end
% 
% attached_eNodeB_pos = BS_pos(eNodeB_order,:);
% angle = mod(atan2(UE_pos(2)-attached_eNodeB_pos(2), UE_pos(1)-attached_eNodeB_pos(1))* 180 / pi, 360);
% 
% for k = 1:length(attached_sectors)
%     az = attached_sectors(k).azimuth;
%     up_az = az + 60;
%     low_az = az - 60;
%     if az == 0
%         low_az = 300;
%         if angle <= up_az || angle >=low_az
%             eNodeB_sector_id = (eNodeB_order-1)*3+k;
%             break;
%         end
%     end
%     if angle <= up_az && angle >=low_az
%         eNodeB_sector_id = (eNodeB_order-1)*3+k;
%         break;
%     end
% end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
