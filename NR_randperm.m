function P = NR_randperm(N,K)
% �Ͱ汾��matlab��֧��randperm(N,K)������ʽ������
% �������
% N������һ��1��N���������У�������
% K��������������ȡǰK��
% �������
% P��һ������˳���˵�1��N���������У�ȡǰK��

P_1 = randperm(N);
try
    P = P_1(1:K);
catch
    % �������K����N�����
    error('K must <= N');
end

end

