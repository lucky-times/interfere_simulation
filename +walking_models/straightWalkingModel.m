classdef straightWalkingModel < walking_models.walkingModel
    % UEֱ������ģ��
    
    properties
        % ����
        direction
        % ��UE_speed*TTI_length���·��
        s
    end
    methods
        % ���캯��
        function obj = straightWalkingModel(s,varargin)
            % ���û��ָ�����������ѡ��һ������
            if length(varargin)<1
                direction = floor(random('unif',0,359));
            else
                direction = varargin{1};
            end
            obj.direction = direction;
            obj.s = s;
        end
        
        % �������ڵ�λ�������һ��TTI��λ��
        function new_pos = move(obj,current_pos)
            mov_vector = obj.s * [ cosd(obj.direction) sind(obj.direction) ];
            new_pos = current_pos + mov_vector;
        end
        
        % �������ڵ�λ�������һ��TTI��λ��
        function old_pos = move_back(obj,current_pos)
            mov_vector = obj.s * [ cosd(obj.direction) sind(obj.direction) ];
            old_pos = current_pos - mov_vector;
        end
        
        % Print some info
        function print(obj)
            fprintf('  direction: %d? %d m/s*TTI\n',obj.direction,obj.s);
        end
    end
end
