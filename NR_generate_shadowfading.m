function networkShadowFadingMap = NR_generate_shadowfading(networkPathlossMap,eNodeB_sites,SYS_config)
%���ڲ�����Ӱ˥��
%���������
%networkPathlossMap��·��ͼ
%eNodeBs��eNodeBʵ��
%SYS_config����������
%
%���������
%networkShadowFadingMap����Ӱ˥��ͼ��
if SYS_config.debug_level>=1
    fprintf('Generating shadow fading\n');
end

networkShadowFadingMap = channel_gain_wrappers.shadowFadingMapClaussen(SYS_config,networkPathlossMap,eNodeB_sites);

if SYS_config.debug_level>=2
    networkShadowFadingMap.print;
end

end