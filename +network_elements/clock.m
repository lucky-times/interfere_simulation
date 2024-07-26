classdef clock < handle
% ʱ���࣬����network element ����һ��ʱ��

   properties
       current_TTI  % TTI��
       TTI_time     % һ��TTI������ʱ��
       time         % ʵ��ʱ��
   end

   methods
       % ���캯��
       function obj = clock(TTI_time)
           obj.current_TTI = 0;
           obj.time        = 0;
           obj.TTI_time    = TTI_time;
       end
       
       % TTI��һ
       function advance_1_TTI(obj)
           obj.current_TTI = obj.current_TTI + 1;
           obj.time = obj.time + obj.TTI_time;
       end
   end
end 
