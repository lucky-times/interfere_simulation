function activate_UE_id = NR_get_activate_UE( UEs )
% �����Ҫ����ͳ�Ƶ�UE��id
% �������
% activate_UE_id����Ҫ����ͳ�Ƶ�UE��id

activate_UE_id = [];
for u_ = 1:length(UEs)
    if ~UEs(u_).deactivate_UE
        if isempty(activate_UE_id)
            activate_UE_id = UEs(u_).id;
        else
            activate_UE_id = [activate_UE_id UEs(u_).id];
        end
    end
end

end

