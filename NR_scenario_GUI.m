function varargout = NR_scenario_GUI(varargin)
% NR_SCENARIO_GUI MATLAB code for NR_scenario_GUI.fig
%      NR_SCENARIO_GUI, by itself, creates a new NR_SCENARIO_GUI or raises the existing
%      singleton*.
%
%      H = NR_SCENARIO_GUI returns the handle to a new NR_SCENARIO_GUI or the handle to
%      the existing singleton*.
%
%      NR_SCENARIO_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NR_SCENARIO_GUI.M with the given input arguments.
%
%      NR_SCENARIO_GUI('Property','Value',...) creates a new NR_SCENARIO_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NR_scenario_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NR_scenario_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NR_scenario_GUI

% Last Modified by GUIDE v2.5 23-Jun-2017 13:51:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NR_scenario_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @NR_scenario_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before NR_scenario_GUI is made visible.
function NR_scenario_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NR_scenario_GUI (see VARARGIN)

% Choose default command line output for NR_scenario_GUI
% �����ʼ������
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.homo,'Value',1); % ��ʼ��Ĭ��Ϊͬ���龰
set(handles.UMA,'Value',1); % ��ʼ��Ĭ��ΪUMA����
set(handles.ue_bf,'Value',1); % ��ʼ��Ĭ��UE����BF
set(handles.bs_bf,'Value',1); % Ĭ��BSʹ��BF
set(handles.issave,'Value',1); % Ĭ�ϱ��沨����������
set(handles.use_cache,'Value',1); % Ĭ��ʹ�ò������ͻ���
set(handles.still,'Value',1); % Ĭ��UE�Ǿ�ֹ��
set(handles.wrap,'Value',1); %Ĭ�ϲ�ʹ��wraparound

set(handles.Rdrop,'Enable','off'); % Ĭ��ͬ��ϵͳΪUMAϵͳ���������޸�UMI������
set(handles.femto,'Enable','off'); % Ĭ��ͬ��ϵͳΪUMAϵͳ���������޸�InH������
set(handles.UMA2UMI,'Enable','off'); % ������ѡ��UMA_to_UMI���칹����
set(handles.UMI2InH,'Enable','off'); % ������ѡ��UMI_to_InH���칹����
set(handles.UMA2InH,'Enable','off'); % ������ѡ��UMA_to_InH���칹����
set(handles.F_2_L,'Enable','off'); % ������ѡ��ǰ�߸��ź���
set(handles.L_2_F,'Enable','off'); % ������ѡ����߸���ǰ��
set(handles.bs_max_gain2,'Enable','off'); % �����������칹���Ԫ����
set(handles.h_bs2,'Enable','off'); % �����������칹��վ�߶�
set(handles.p_bs2,'Enable','off'); % �����������칹��վ���书��
set(handles.speed,'Enable','off'); % �����������ٶȣ���Ϊ�û��ƶ�û��ѡ��
set(handles.TTI,'Enable','off'); % �������޸ķ���ʱ������Ϊ�û��ƶ�û��ѡ��
set(handles.div_angle,'Enable','off'); % ����������ƫ�ƽǶ�
set(handles.div_dist,'Enable','off'); % ����������ƫ�ƾ���


set(handles.bs_array_out,'String','[8,16]'); % Ĭ�ϻ�������������Ϊ8*16
set(handles.bs_array_in,'String','[4,8]'); % Ĭ����������������Ϊ4*8
set(handles.ue_max_gain,'String','5'); % Ĭ���û��������Ԫ����Ϊ5dB
set(handles.bs_max_gain1,'String','8'); % Ĭ��victimϵͳ����ͬ����aggressorϵͳ����վ�������Ԫ����Ϊ8dB
set(handles.bs_phi,'String','0'); % Ĭ�ϻ�վ�ನ����ˮƽ������Ͻ�Ϊ0
set(handles.bs_theta,'String','0'); % Ĭ�ϻ�վ�ನ����ֱ������Ͻ�Ϊ0
set(handles.ue_phi,'String','0'); % Ĭ��UE�ನ��ˮƽ������Ͻ�Ϊ0
set(handles.ue_theta,'String','0'); % Ĭ��UE�ನ����ֱ������Ͻ�Ϊ0
set(handles.acir_lower,'String','5'); % Ĭ��ACIR�о��½�Ϊ5dB
set(handles.acir_upper,'String','40'); % Ĭ��ACIR�о��Ͻ�Ϊ40dB
set(handles.acir_step,'String','5'); % Ĭ��ACIR�о�����ֵΪ5dB
set(handles.acir,'String','15'); % Ĭ��ACIRУ׼ֵΪ15dB
set(handles.slot,'String','5'); % Ĭ�Ͽ�����Ϊ10
set(handles.f1,'String','30'); % Ĭ��victimϵͳƵ��Ϊ30GHz
set(handles.f2,'String','30.2'); % Ĭ��aggressorϵͳƵ��Ϊ30.2GHz
set(handles.band,'String','200'); % Ĭ�ϴ���Ϊ200MHz
set(handles.p_ue,'String','23'); % Ĭ��UE�ķ��书��Ϊ23dBm
set(handles.p_bs1,'String','43'); % Ĭ��victimϵͳ����ͬ����aggressorϵͳ����վ���书��Ϊ43dBm
set(handles.nf_ue,'String','11'); % Ĭ��UE������ϵ��Ϊ11dB
set(handles.nf_bs,'String','11'); % Ĭ��BS������ϵ��Ϊ11dB
set(handles.h_bs1,'String','25'); % Ĭ��victimϵͳ����ͬ����aggressorϵͳ����վ�߶�Ϊ25m
set(handles.shadow_r,'String','0.5'); % Ĭ����Ӱ˥�����ϵ��Ϊ0.5
set(handles.res,'String','5'); % Ĭ�Ϸֱ���Ϊ5
set(handles.isd,'String','200'); % Ĭ��ϵͳISDΪ200m
set(handles.speed,'String','0'); %Ĭ���û��ٶ�Ϊ0
set(handles.TTI,'String','1'); % Ĭ�Ϸ���ʱ��Ϊ1��TTI
set(handles.UEs_per_slot,'String','1'); % Ĭ��ÿ�ο��շ���1���û�


% UIWAIT makes NR_scenario_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NR_scenario_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in structure.
function structure_Callback(hObject, eventdata, handles)
% hObject    handle to structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns structure contents as cell array
%        contents{get(hObject,'Value')} returns selected item from structure


% --- Executes during object creation, after setting all properties.
function structure_CreateFcn(hObject, eventdata, handles)
% hObject    handle to structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hete.
function hete_Callback(hObject, eventdata, handles)
% hObject    handle to hete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% �칹��ť������

set(handles.homo,'Value',0); %����ͬ��ѡ��
%���þ���ͬ��������ѡ��
set(handles.RMA,'Value',0);  
set(handles.UMA,'Value',0);
set(handles.InH,'Value',0);
set(handles.UMI,'Value',0);
set(handles.F_2_L,'Enable','on');%����ѡ���������
set(handles.L_2_F,'Enable','on');%����ѡ�������
set(handles.F_2_L,'Value',1);%Ĭ���칹������ǰ�߸��ź��ߣ�����
set(handles.L_2_F,'Value',0);
%Ĭ���칹����ΪUMA_to_UMI
set(handles.UMA2UMI,'Value',1);
set(handles.loc,'Value',2); %�칹����Ĭ��ѡ�񲻹�ַ
%��������ͬ��������ѡ��
set(handles.RMA,'Enable','off');
set(handles.UMA,'Enable','off');
set(handles.InH,'Enable','off');
set(handles.UMI,'Enable','off');
set(handles.Rdrop,'Enable','on'); % Ĭ���칹ϵͳΪUMa_to_UMi.�ʽ���UMI���ⳡ��ѡ��֧��UMA��manhattan��UMA��Rdrop��������
set(handles.Rdrop,'Value',1);
set(handles.femto,'Enable','off'); % ����InH����ѡ��ڴ˳���������
set(handles.wrap,'Enable','off'); %����wraparoundѡ��
% �����칹��������
set(handles.UMA2UMI,'Enable','on'); 
set(handles.UMA2InH,'Enable','on');
set(handles.UMI2InH,'Enable','on');
set(handles.bs_max_gain2,'Enable','on'); %�����칹ϵͳ��վ�������Ԫ���������
set(handles.h_bs2,'Enable','on'); %�����칹ϵͳ��վ�߶ȵ����� 
set(handles.p_bs2,'Enable','on'); %�����칹ϵͳ��վ���书�ʵ�����
set(handles.loc,'Enable','off'); %������ַ��ѡ�񣬽�ֹѡ��Ϊ��ַѡ��
% set(handles.together,'Enable','off'); %����ͬ��/�첽�������͡����칹ϵͳ��֧���첽���ŵ����
set(handles.div_angle,'Enable','off'); %����С����ƫ��󳡾����ĵĽǶ�
set(handles.div_dist,'Enable','off'); %����С����ƫ��󳡾����ĵľ���
%Ĭ��С�����ڴ󳡾�������
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.sysnum,'Enable','off'); %����ϵͳ��ѡ����칹ϵͳ��ϵͳ��Ϊ2
set(handles.isd,'Enable','on'); %����ISD���ã�������������victimϵͳ��ISD
%����ACIR��ص����á�����ϵͳѡ��ʱ�ὫACIR���ֵ��ղ�������
set(handles.acir_upper,'Enable','on'); 
set(handles.acir_lower,'Enable','on');
set(handles.acir_step,'Enable','on');
set(handles.acir,'Enable','on');
%����ϵͳ����ѡ�񡪡�ֻ��UMA/RMA֧��LTE��NR�Ĺ������
set(handles.v_type,'Enable','off');
set(handles.a_type,'Enable','on');


set(handles.bs_max_gain2,'String','8'); %���������칹��վ�������Ԫ����Ϊ8dB
set(handles.h_bs2,'String','25'); %���������칹ϵͳ��վ�߶�Ϊ25m
set(handles.p_bs2,'String','43'); %���������칹ϵͳ��վ���书��Ϊ43dBm
set(handles.p_bs1,'String','33'); %��������victimϵͳ��վ���书��Ϊ33dBm
set(handles.h_bs1,'String','10'); %��������victimϵͳ��վ�߶�Ϊ10m
set(handles.slot,'String','5'); %��������ϵͳ������Ϊ5
 %������ϵͳƵ��Ϊ30G��30.2GHz����RMa�����ᴥ��Ƶ�ʵĸ���
set(handles.f1,'String','30'); 
set(handles.f2,'String','30.2');
%��������ϵͳ����ͬ������ģʽ
set(handles.together,'Value',1);
set(handles.sysnum,'Value',1); %��������ϵͳ��Ϊ2
set(handles.wrap,'Value',1); %�������ò�����wraparound����ֻ��UMA/RMA����֧��wraparound
set(handles.isd,'String','200'); %��������ϵͳISD��200
%��������ACIR��ز��������á�������ӵ�ϵͳ״̬�������칹״̬������Ҫ����ACIR��ز�������ֵ
set(handles.acir,'String','15');
set(handles.acir_step,'String','5');
set(handles.acir_lower,'String','5');
set(handles.acir_upper,'String','40');
set(handles.handover,'Value',1); % �������ò���3dB handover
set(handles.shadow_r,'String','0.5'); % ����������Ӱ˥�����ϵ��Ϊ0.5
%������ϵͳ������ΪNR
set(handles.a_type,'Value',1);
set(handles.v_type,'Value',1);
% Hint: get(hObject,'Value') returns toggle state of hete


% --- Executes on button press in homo.
function homo_Callback(hObject, eventdata, handles)
% hObject    handle to homo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ���ͬ��ѡ��Ĵ����¼�
% Ĭ��ΪUMAϵͳ�����������ֳ���������
set(handles.Rdrop,'Enable','off');
set(handles.femto,'Enable','off');
%�����칹ϵͳ���峡����ѡ�񲢿���ֵ
set(handles.UMA2UMI,'Value',0);
set(handles.UMI2InH,'Value',0);
set(handles.UMA2InH,'Value',0);
set(handles.UMA2UMI,'Enable','off');
set(handles.UMI2InH,'Enable','off');
set(handles.UMA2InH,'Enable','off');
set(handles.F_2_L,'Enable','off');%������ѡ���������
set(handles.L_2_F,'Enable','off');%������ѡ�������
%�칹���ŷ������λ
set(handles.F_2_L,'Value',0);
set(handles.L_2_F,'Value',0);
%��ͬ��ϵͳ���峡����ѡ��
set(handles.UMA,'Enable','on');
set(handles.RMA,'Enable','on');
set(handles.UMI,'Enable','on');
set(handles.InH,'Enable','on');
%�����칹ϵͳר�в��������Ԫ���棬���߸߶ȣ����书�ʣ�
set(handles.bs_max_gain2,'Enable','off');
set(handles.h_bs2,'Enable','off');
set(handles.p_bs2,'Enable','off');
%�����첽����ѡ��
set(handles.together,'Enable','on');
set(handles.sysnum,'Enable','on'); %������ϵͳѡ��
%�����칹ϵͳĬ��ѡ��UMA����
set(handles.loc,'Enable','on'); %��������ַѡ��
set(handles.wrap,'Enable','on'); %����wraparound��ѡ��


set(handles.bs_array_out,'String','[8,16]');
set(handles.ue_max_gain,'String','5');
set(handles.bs_max_gain1,'String','8'); %���û�վ���������Ԫ����
set(handles.slot,'String','5');% ���ÿ�����
set(handles.h_bs1,'String','25'); % ���û�վ�����߸߶�
set(handles.res,'String','5'); %���÷ֱ���
set(handles.isd,'String','200'); %����ISD
%����칹ϵͳר�в���
set(handles.h_bs2,'String','');  
set(handles.p_bs2,'String','');
set(handles.bs_max_gain2,'String','');

set(handles.together,'Value',1); %Ĭ��Ϊͬ������
set(handles.hete,'Value',0); %����ͬ��ѡ���ͬ���칹ֻ��ѡ��һ��
set(handles.sysnum,'Value',1); %����ϵͳ��Ϊ˫ϵͳ
set(handles.UMA,'Value',1); %Ĭ��ͬ������ΪUMA
set(handles.loc,'Value',1); %���ù�ַѡ��
%�����ϵͳƫ�ýǶ���ƫ��ֵ
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
set(handles.wrap,'Value',1); %Ĭ�ϲ�����wraparound
%����ƫ�ýǶ�����롪��Ĭ��ͬ��ϵͳΪUMA
set(handles.v_type,'Enable','on'); 
set(handles.a_type,'Enable','on');

% Hint: get(hObject,'Value') returns toggle state of homo


% --- Executes on selection change in sysnum.
function sysnum_Callback(hObject, eventdata, handles)
% hObject    handle to sysnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%������/˫ϵͳ������

n_sysnum = get(hObject,'value'); %�õ�ѡ������ݣ�1Ϊ˫ϵͳ��2Ϊ��ϵͳ
if n_sysnum==1 %˫ϵͳ
    set(handles.loc,'Value',1); %Ĭ�Ϲ�ַ�龰
    if get(handles.UMA,'Value') + get(handles.RMA,'Value')>0 %�������UMA��RMA
        set(handles.loc,'Enable','on'); %��������ַѡ�� ����ֻ��UMA��RMA֧����������ƫ��״̬
    else
        set(handles.loc,'Enable','off'); %��������ַѡ��
    end
    %����ACIR������ò��ָ�ACIR�����Ĭ������
    set(handles.acir_lower,'String','5');
    set(handles.acir_upper,'String','40');
    set(handles.acir,'String','15');
    set(handles.acir_step,'String','5');
    set(handles.acir_upper,'Enable','on');
    set(handles.acir_lower,'Enable','on');
    set(handles.acir_step,'Enable','on');
    set(handles.acir,'Enable','on');
    %��������ģʽѡ��Ĭ��Ϊͬ������
    set(handles.together,'Value',1);
    set(handles.together,'Enable','on');
else %��ϵͳ
    %���������ù�ַѡ��
    set(handles.loc,'Value',0);
    set(handles.loc,'Enable','off');
    %����������ƫ�ýǶ���ƫ�þ���
    set(handles.div_angle,'Enable','off');
    set(handles.div_dist,'Enable','off');
    set(handles.div_angle,'String','');
    set(handles.div_dist,'String','');
    %����������ACIR���������
    set(handles.acir_lower,'String','');
    set(handles.acir_upper,'String','');
    set(handles.acir,'String','');
    set(handles.acir_step,'String','');
    set(handles.acir_upper,'Enable','off');
    set(handles.acir_lower,'Enable','off');
    set(handles.acir_step,'Enable','off');
    set(handles.acir,'Enable','off');
    %��������ģʽ��ѡ�񣬲�����Ϊͬ������
    set(handles.together,'Enable','off');
    set(handles.together,'Value',0);
end
   
% Hints: contents = cellstr(get(hObject,'String')) returns sysnum contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sysnum


% --- Executes during object creation, after setting all properties.
function sysnum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sysnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in OK.
function OK_Callback(hObject, eventdata, handles)

% hObject    handle to OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ���ȷ����ť�Ĵ����¼�
% ���ڲ�����z_test����������еı����������GUI�ľ������ֻ�ܲ���clc��ղ������ݡ���ͼ�ν������г�����ű�����������
% ������clear��������GUI
% �ô����¼���Ҫ�ǽ��в����Ĵ��ݺ�������������GUI�ĵ���
% �Ӹú����п��Ժ�����ؿ����������������Ķ�Ӧ��ϵ
clc;
SYS_config = NR_load_params;%����Ĭ�Ϸ������
%% ����Խ����ı�����Ĭ�����ã�����ͼ�ν���������
SYS_config.shadow_fading_type         = 'claussen';
SYS_config.compact_results_file       = true;
SYS_config.TTI_length = 1;% Length of a TTI (subframe), in seconds.
SYS_config.default_shown_GUI_cells = [];
SYS_config.PCe_param = 0; %���صľ�ȷ��
SYS_config.interference_type = 2;% 0ͬƵ���ţ�1��Ƶ���ţ�2���У�ֻ��˫ϵͳ����Ч
SYS_config.attatch_mode=3;  % 1��ʾ·���BS�˵�Ԫ���棻2��ʾ·���UE�˵�Ԫ���棻3��ʾ·���BS��UE�˵�Ԫ���档����ȷ��UE����С�������㣩
SYS_config.Gama = 1; %�����е�gamma����
SYS_config.antenna.antenna_gain_pattern = 'NRAntennaBeamforming'; %�������͵����ͣ������ڴ˴�������չ�����Ĳ�������ͼ����
SYS_config.AntennaPattern3d = true;
SYS_config.beam_loss=0; %�������
SYS_config.cable_loss=0; %�����������

SYS_config.ISD = str2double(get(handles.isd,'String'));
%% ר�ò���
if get(handles.UMA,'Value')==1 %UMA
    SYS_config.scene_type = 'UMA';
    SYS_config.isNewUMA = true;
    SYS_config.sector_azimuths = 0:360/3:359;
    SYS_config.shadow_fading_sd_LOS = 4; %LOS����Ӱ˥���׼��
    SYS_config.shadow_fading_sd_NLOS = 6; %NLOS����Ӱ˥���׼��
    SYS_config.default_shown_GUI_cells = 1:57;
    SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_macro'; %��һ��ϵͳ�Ļ���
    SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro'; %�ڶ���ϵͳ�Ļ���
end
if get(handles.UMI,'Value')==1 %UMI
    SYS_config.scene_type = 'UMI';
    SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
    SYS_config.sector_azimuths = 0;% ��ʼ�������Ƕ�
    SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_micro';
    SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_micro';
    SYS_config.shadow_fading_sd_LOS = 4;
    SYS_config.shadow_fading_sd_NLOS = 7.82;
end
if get(handles.RMA,'Value')==1 %RMA
    SYS_config.scene_type = 'RMa';
    SYS_config.sector_azimuths = 0:360/3:359;
    SYS_config.macroscopic_pathloss_model_settings.environment = 'rural_macro';
    SYS_config.macroscopic_pathloss_model_settings2.environment = 'rural_macro';
    SYS_config.shadow_fading_sd_LOS = 4;
    SYS_config.shadow_fading_sd_LOS2 = 6; %900��RMA������LOS״̬������LOS��
    SYS_config.shadow_fading_sd_NLOS = 8;
    SYS_config.default_shown_GUI_cells = 1:57;
end
if get(handles.InH,'Value')==1 %InH
    SYS_config.scene_type = 'InH';
    SYS_config.sector_azimuths = 0;% ��ʼ�������Ƕ�
    SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
    SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
    SYS_config.shadow_fading_sd_LOS = 3;
    SYS_config.shadow_fading_sd_NLOS = 8.03;
    SYS_config.UE_height=1;
end
if get(handles.UMA2UMI,'Value')==1 %UMa_to_UMi
    SYS_config.scene_type = 'UMa_to_UMi';
    SYS_config.isS2F = get(handles.F_2_L,'Value');
    SYS_config.isNewUMA = true;
    SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
    SYS_config.sector_azimuths = 0;% ��ʼ�������Ƕ�
    if SYS_config.isS2F %�����ǰ�߸��ź��ߣ���Ϊ����ϵͳ����UMA������������ϵͳ����INH����
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_micro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro';
    else %����Ϊ����ϵͳ����indoor������Ϊ������ϵͳ����UMA����
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_macro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_micro';
    end
    SYS_config.shadow_fading_sd_LOS = 4;
    SYS_config.shadow_fading_sd_NLOS = 7.82;
    SYS_config.shadow_fading_sd_LOS_het = 4;
    SYS_config.shadow_fading_sd_NLOS_het = 6;
    SYS_config.sector_azimuths2 = 0:360/3:359;
   
end
if get(handles.UMA2InH,'Value')==1 %UMa_to_InH
    SYS_config.scene_type = 'UMa_to_InH';
    SYS_config.isS2F = get(handles.F_2_L,'Value');
    SYS_config.sector_azimuths = 0;% ��ʼ�������Ƕ�
    if SYS_config.isS2F %�����ǰ�߸��ź��ߣ���Ϊ����ϵͳ����UMA������������ϵͳ����INH����
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_macro';
    else %����Ϊ����ϵͳ����indoor������Ϊ������ϵͳ����UMA����
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_macro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
    end
    SYS_config.shadow_fading_sd_LOS = 3;
    SYS_config.shadow_fading_sd_NLOS = 8.03;
    SYS_config.sector_azimuths2 = 0:360/3:359;
    SYS_config.shadow_fading_sd_LOS_het = 4;
    SYS_config.shadow_fading_sd_NLOS_het = 6;
    SYS_config.UE_height=1.5;
end
if get(handles.UMI2InH,'Value')==1 %UMi_to_InH
    SYS_config.scene_type = 'UMi_to_InH';
    SYS_config.isS2F = get(handles.F_2_L,'Value');
    SYS_config.sector_azimuths = 0;% ��ʼ�������Ƕ�
    if SYS_config.isS2F %�����ǰ�߸��ź��ߣ���Ϊ����ϵͳ����UMI������������ϵͳ����INH����
        SYS_config.macroscopic_pathloss_model_settings.environment = 'indoor';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'urban_micro';
    else %����Ϊ����ϵͳ����indoor������Ϊ������ϵͳ����UMI����
        SYS_config.macroscopic_pathloss_model_settings.environment = 'urban_micro';
        SYS_config.macroscopic_pathloss_model_settings2.environment = 'indoor';
    end
    SYS_config.shadow_fading_sd_LOS = 3;
    SYS_config.shadow_fading_sd_NLOS = 8.03;
    SYS_config.UMi_r = (SYS_config.ISD/3)/2/2*(3^0.5);
    SYS_config.sector_azimuths2 = 0;
    SYS_config.shadow_fading_sd_LOS_het = 4;
    SYS_config.shadow_fading_sd_NLOS_het = 7.82;
    SYS_config.UE_height=1.5;
end
if get(handles.Rdrop,'Value')==1
    SYS_config.isManhattan = false;
else
    SYS_config.isManhattan = true;
end
if get(handles.femto,'Value')==1
    SYS_config.isFemto = false;
else
    SYS_config.isFemto = true;
end

SYS_config.simulation_time_tti = str2double(get(handles.TTI,'String'));
SYS_config.keep_UEs_still  = get(handles.still,'Value')==1;
SYS_config.UE_speed = str2double(get(handles.TTI,'String'))/3.6; %���ٶ�ֵת��Ϊm/s

SYS_config.beam_accuracy_phi_ue = str2double(get(handles.ue_phi,'String'));
SYS_config.beam_accuracy_theta_ue = str2double(get(handles.ue_theta,'String'));
SYS_config.beam_accuracy_phi_bs = str2double(get(handles.bs_phi,'String'));
SYS_config.beam_accuracy_theta_bs = str2double(get(handles.bs_theta,'String'));

if get(handles.distribution,'Value')==1
    SYS_config.beam_accurancy_type = 'Constant';
else
    SYS_config.beam_accurancy_type = 'Unique';
end

SYS_config.asynchronization = get(handles.together,'Value')==2;


SYS_config.frequency = str2double(get(handles.f1,'String')) * 1e9; %��λת������GHzת��ΪHz
SYS_config.frequency2 = str2double(get(handles.f2,'String')) * 1e9; %��λ���껻����GHzת��ΪHz
SYS_config.bandwidth = str2double(get(handles.band,'String')) * 1e6; %��λת������MHzת��ΪHz
SYS_config.photo = get(handles.isphoto,'Value');
SYS_config.ACIR_dB = str2double(get(handles.acir,'String'));
SYS_config.isACIR_loop = true;
SYS_config.ACIR_lower = str2double(get(handles.acir_lower,'String'));
SYS_config.ACIR_upper = str2double(get(handles.acir_upper,'String'));
SYS_config.isDouble = get(handles.sysnum,'Value')==1;
if get(handles.ue_bf,'value')==1 && get(handles.bs_bf,'value')==1
    SYS_config.antenna_mode = 2;
elseif get(handles.bs_bf,'value')==1
    SYS_config.antenna_mode = 1;
else
    SYS_config.antenna_mode = 0;
end
SYS_config.use_cache = get(handles.use_cache,'value')==1;
SYS_config.isSave = get(handles.issave,'value')==1;
SYS_config.antenna_element = str2num(get(handles.bs_array_out,'String'));
SYS_config.antenna_element_InH = str2num(get(handles.bs_array_in,'String'));
SYS_config.UE_tx_power = 10^(0.1*str2double(get(handles.p_ue,'String')))/1000;
SYS_config.UE_max_antenna_gain = str2double(get(handles.ue_max_gain,'String'));
SYS_config.UE_receiver_noise_figure = str2double(get(handles.nf_ue,'String'));
SYS_config.BS_receiver_noise_figure = str2double(get(handles.nf_bs,'String'));
if get(handles.v_type,'Value')==1
    SYS_config.macroscopic_pathloss_model = 'TS38901'; %901��ӦNRϵͳ��·��ģ��
else
    SYS_config.macroscopic_pathloss_model = 'TS36942'; %942��ӦLTEϵͳ��·��ģ��
end
if get(handles.a_type,'Value')==1
    SYS_config.macroscopic_pathloss_model2 = 'TS38901';
else
    SYS_config.macroscopic_pathloss_model2 = 'TS36942';
end
SYS_config.isWraparound = get(handles.wrap,'Value')==2;
SYS_config.ISD = str2double(get(handles.isd,'String'));
if get(handles.loc,'Value')==1
    SYS_config.shift_mode = 0;
else
    SYS_config.shift_mode = 2;
%     if get(handles.UMA,'Value') + get(handles.RMA,'Value')==0 
%         SYS_config.shift_mode = 1;
%     else
%         SYS_config.shift_mode = 2;
%     end
end
SYS_config.angle_between_systems = str2double(get(handles.div_angle,'String'));
SYS_config.ISD_between_systems = str2double(get(handles.div_dist,'String'));
SYS_config.hand_over = get(handles.handover,'Value')==1;
SYS_config.map_resolution = str2double(get(handles.res,'String'));
SYS_config.shadow_fading_map_resolution = str2double(get(handles.res,'String'));
SYS_config.n_snapshot = str2double(get(handles.slot,'String'));
SYS_config.n_snapshot2 = str2double(get(handles.slot,'String'));
SYS_config.n_UE_served_per_BS = str2double(get(handles.UEs_per_slot,'String'));
SYS_config.r_eNodeBs = str2double(get(handles.shadow_r,'String'));
SYS_config.site_height = str2double(get(handles.h_bs1,'String'));
SYS_config.eNodeB_tx_power = 10^(0.1*str2double(get(handles.p_bs1,'String')))/1000; %��dBmת��ΪW
SYS_config.antenna.max_antenna_gain = str2double(get(handles.bs_max_gain1,'String')); %���ַ�ת��Ϊ��ֵ
SYS_config.eNodeB_tx_power2 = 10^(0.1*str2double(get(handles.p_bs2,'String')))/1000; %��dBmת��ΪW
SYS_config.site_height2 = str2double(get(handles.h_bs2,'String'));
SYS_config.antenna.max_antenna_gain2 = str2double(get(handles.bs_max_gain2,'String'));

 
 output_results_file = NR_sim_main(SYS_config); %����������
 simulation_data                   = load(output_results_file); %��������
 
 GUI_handles.positions_GUI         = LTE_GUI_show_UEs_and_cells(simulation_data); %��������ͼ
% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
% hObject    handle to cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ȡ����ť������ʱ��
clear; %������еı���
close(gcf); %�ر�GUI

% --- Executes on button press in UMA2InH.
function UMA2InH_Callback(hObject, eventdata, handles)
% hObject    handle to UMA2InH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����UMA����InH����������
% InH���Ϳ�ѡ�����漰UMI�������������Ͳ���ѡ
set(handles.Rdrop,'Enable','off');
set(handles.Rdrop,'Value',1);
set(handles.femto,'Enable','on');
% ���������칹������ѡ��״̬
set(handles.UMA2UMI,'Value',0);
set(handles.UMI2InH,'Value',0);
%������ϵͳ��վ����
set(handles.p_bs1,'String',23);
set(handles.h_bs1,'String',3);
set(handles.p_bs2,'String',43);
set(handles.h_bs2,'String',25);
set(handles.bs_max_gain1,'String',5);
set(handles.res,'String',2); %���÷ֱ���
set(handles.slot,'String',5); %���ÿ�����
%������Ӱ˥�����ϵ��
set(handles.shadow_r,'String','0.5');
set(handles.v_type,'Value',1);%Ĭ�ϱ�����ϵͳΪNR
set(handles.a_type,'Value',1);%Ĭ�ϸ���ϵͳΪNR
if get(handles.F_2_L,'Value') %�������ϵͳ��UMA
    set(handles.v_type,'Enable','off'); %������ϵͳһ��ΪNR�Ҳ��ɱ�
    set(handles.a_type,'Enable','on'); %����ϵͳ��ѡLTE
else  %���������ϵͳ��UMA
    set(handles.v_type,'Enable','on');%������ϵͳ��ѡLTE
    set(handles.a_type,'Enable','off');%����ϵͳһ��ΪNR
end
%����ƫ�÷���Ǻ�ƫ�þ���
set(handles.div_angle,'Enable','on');
set(handles.div_dist,'Enable','on');
set(handles.div_angle,'String','0');
set(handles.div_dist,'String','0');
% Hint: get(hObject,'Value') returns toggle state of UMA2InH


% --- Executes on button press in UMI2InH.
function UMI2InH_Callback(hObject, eventdata, handles)
% hObject    handle to UMI2InH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����UMI����InH����������
% ���ڲ�֧��manhattan��InH���칹�������ʸó�����InH���Ϳ�ѡ��UMI��������ѡ
set(handles.Rdrop,'Value',1); %����Ϊrandom drop�����о���InH�����Ĺ���
set(handles.Rdrop,'Enable','off');
set(handles.femto,'Enable','on');
% ���������칹������ѡ��״̬
set(handles.UMA2UMI,'Value',0);
set(handles.UMA2InH,'Value',0);
%������ϵͳ��վ����
set(handles.p_bs1,'String',23);
set(handles.h_bs1,'String',3);
set(handles.p_bs2,'String',33);
set(handles.h_bs2,'String',10);
set(handles.bs_max_gain1,'String',5);
set(handles.res,'String',2); %���÷ֱ���
set(handles.slot,'String',5); %���ÿ�����
set(handles.shadow_r,'String','0.5'); %������Ӱ˥�����ϵ��
%��ϵͳ��ΪNR�Ҳ���ѡLTE
set(handles.v_type,'Enable','off');
set(handles.a_type,'Enable','off');
set(handles.v_type,'Value',1);
set(handles.a_type,'Value',1);
% Hint: get(hObject,'Value') returns toggle state of UMI2InH
%����ƫ�÷���Ǻ�ƫ�þ���
set(handles.div_angle,'Enable','on');
set(handles.div_dist,'Enable','on');
set(handles.div_angle,'String','0');
set(handles.div_dist,'String','0');

% --- Executes on button press in UMA2UMI.
function UMA2UMI_Callback(hObject, eventdata, handles)
% hObject    handle to UMA2UMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����UMA����UMI���������á�
% ƽ̨֧��UMA��manhattan��UMA��Rdrop����������칹��������UMI���Ϳ�ѡ�����漰InH������������ѡ��InH����
set(handles.Rdrop,'Enable','on');
set(handles.Rdrop,'Value',1);
set(handles.femto,'Enable','off');
% ���������칹������ѡ��״̬
set(handles.UMA2InH,'Value',0);
set(handles.UMI2InH,'Value',0);
%������ϵͳ��վ����
set(handles.p_bs1,'String',33);
set(handles.h_bs1,'String',10);
set(handles.p_bs2,'String',43);
set(handles.h_bs2,'String',25);
set(handles.bs_max_gain1,'String',8);
set(handles.res,'String',5); %���÷ֱ���
set(handles.slot,'String',5);%���ÿ�����
 %������Ӱ˥�����ϵ��
set(handles.shadow_r,'String','0.5');
set(handles.v_type,'Value',1);
set(handles.a_type,'Value',1);
if get(handles.F_2_L,'Value') %�������ϵͳ��UMA
    set(handles.v_type,'Enable','off'); %������ϵͳһ��ΪNR�Ҳ��ɱ�
    set(handles.a_type,'Enable','on'); %����ϵͳ��ѡLTE
else  %���������ϵͳ��UMA
    set(handles.v_type,'Enable','on');%������ϵͳ��ѡLTE
    set(handles.a_type,'Enable','off');%����ϵͳһ��ΪNR
end
%����ƫ�÷���Ǻ�ƫ�þ���
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
    
% Hint: get(hObject,'Value') returns toggle state of UMA2UMI


% --- Executes on button press in InH.
function InH_Callback(hObject, eventdata, handles)
% hObject    handle to InH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ѡ��InH�����µĴ����¼�
% ��������ͬ��������ѡ��״̬
set(handles.RMA,'Value',0);
set(handles.UMA,'Value',0);
set(handles.UMI,'Value',0);
%����wraparoundѡ�����Ϊ������wraparound
set(handles.wrap,'Value',1);
set(handles.wrap,'Enable','off');
set(handles.Rdrop,'Enable','off');%����UMI����������
set(handles.femto,'Enable','on');%����InH����������
set(handles.slot,'String','5'); %���ÿ�����
%���û�վ��ز���
set(handles.p_bs1,'String','23');
set(handles.h_bs1,'String','3');
set(handles.bs_max_gain1,'String','5');
%��ղ�����ƫ��״̬������
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
% set(handles.loc,'Value',1);
% set(handles.loc,'Enable','on');
set(handles.res,'String','0.5'); %���÷ֱ���
%����Ƶ�ʵ�����
set(handles.f1,'String','30');
set(handles.f2,'String','30.2');
%����������ISD���á���InH�����������޸�ISD
set(handles.isd,'String','20');
set(handles.isd,'Enable','off');
set(handles.handover,'Value',2); %����handover״̬Ϊ������3dB handover
set(handles.shadow_r,'String','0'); %������Ӱ˥������ϵ��Ϊ0
if get(handles.sysnum,'Value')==1 %�����˫ϵͳ��رչ�ַ��ѡ������Ĭ��ֵΪ��ַ����InH��֧��ϵͳƫ��
    set(handles.loc,'Value',1);
    set(handles.loc,'Enable','off');
end
% set(handles.sysnum,'Value',1);
% ����������ϵͳ����ѡ��
set(handles.v_type,'Enable','off');
set(handles.a_type,'Enable','off');
set(handles.v_type,'Value',1);
set(handles.a_type,'Value',1);

% Hint: get(hObject,'Value') returns toggle state of InH


% --- Executes on button press in UMI.
function UMI_Callback(hObject, eventdata, handles)
% hObject    handle to UMI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ѡ��UMi�����µĴ����¼�
% ��������ͬ��������ѡ��״̬
set(handles.RMA,'Value',0);
set(handles.UMA,'Value',0);
set(handles.InH,'Value',0);
%����wraparoundѡ�����Ϊ������wraparound
set(handles.wrap,'Value',1);
set(handles.wrap,'Enable','off');
%����UMI�������Ͳ�Ĭ����ΪRandom drop����
set(handles.Rdrop,'Enable','on');
set(handles.femto,'Enable','off');% ����InH����ѡ��
set(handles.Rdrop,'Value',1);
%���û�վ��ز���
set(handles.p_bs1,'String','33');
set(handles.slot,'String','5');
set(handles.h_bs1,'String','10');
set(handles.bs_max_gain1,'String','8');
%��ղ�����ƫ��״̬������
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
% set(handles.loc,'Value',1);
% set(handles.loc,'Enable','on');
set(handles.res,'String','5'); %���÷ֱ���
%����Ƶ�ʵ�����
set(handles.f1,'String','30');
set(handles.f2,'String','30.2');
%����������ISD���á���InH�����������޸�ISD
set(handles.isd,'String','200');
set(handles.isd,'Enable','on');
set(handles.handover,'Value',1); %����handover״̬Ϊ����3dB handover
if get(handles.sysnum,'Value')==1 %�����˫ϵͳ��رչ�ַ��ѡ������Ĭ��ֵΪ����ַ����UMIΪ����ַ����
    set(handles.loc,'Value',2);
    set(handles.loc,'Enable','off');
end
set(handles.shadow_r,'String','0.5');%������Ӱ˥������ϵ��Ϊ0.5
% ����������ϵͳ����ѡ��
set(handles.v_type,'Enable','off');
set(handles.a_type,'Enable','off');
set(handles.v_type,'Value',1);
set(handles.a_type,'Value',1);
% set(handles.sysnum,'Value',1);
% Hint: get(hObject,'Value') returns toggle state of UMI


% --- Executes on selection change in loc.
function loc_Callback(hObject, eventdata, handles)
% hObject    handle to loc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ��ѡ��ַ״̬��ʱ�������¼� 1Ϊ��ַ 2Ϊ����ַ��ֻ��ͬ��ʱ�����޸ģ�
loc_type = get(hObject,'value');
if loc_type==2 %����ַ
    if get(handles.UMA,'Value') + get(handles.RMA,'Value')>0 %RMA����UMA�����������ƫ��״̬��Ĭ��ֵ
        set(handles.div_angle,'Enable','on');
        set(handles.div_dist,'Enable','on');
        set(handles.div_angle,'String','45');
        set(handles.div_dist,'String','50');
    end
else %��ַ�����������ƫ��״̬
    set(handles.div_angle,'Enable','off');
    set(handles.div_dist,'Enable','off');
    set(handles.div_angle,'String','');
    set(handles.div_dist,'String','');
end

% Hints: contents = cellstr(get(hObject,'String')) returns loc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from loc


% --- Executes during object creation, after setting all properties.
function loc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Rdrop.
function Rdrop_Callback(hObject, eventdata, handles)
% hObject    handle to Rdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ��UMI�������͵Ĵ����¼�
if get(handles.Rdrop,'Value')==1 %���ѡ��random drop
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
else %���ѡ��manhattan
    if get(handles.UMA2UMI,'Value')
        set(handles.div_angle,'Enable','on');
        set(handles.div_dist,'Enable','on');
        set(handles.div_angle,'String','0');
        set(handles.div_dist,'String','0');
    end
end
% Hints: contents = cellstr(get(hObject,'String')) returns Rdrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Rdrop


% --- Executes during object creation, after setting all properties.
function Rdrop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function femto_Callback(hObject, eventdata, handles)
% hObject    handle to Rdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ��InH�������͵Ĵ����¼�
% ֻ���ظ�ֵ����
% if get(handles.femto,'Value')==1 %���ѡ��Home femto
% set(handles.div_angle,'Enable','off');
% set(handles.div_dist,'Enable','off');
% set(handles.div_angle,'String','');
% set(handles.div_dist,'String','');
% else %���ѡ��manhattan
%     if get(handles.UMA2UMI,'Value')
%         set(handles.div_angle,'Enable','on');
%         set(handles.div_dist,'Enable','on');
%         set(handles.div_angle,'String','0');
%         set(handles.div_dist,'String','0');
%     end
% end
% Hints: contents = cellstr(get(hObject,'String')) returns Rdrop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Rdrop


% --- Executes during object creation, after setting all properties.
function femto_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Rdrop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in RMA.
function RMA_Callback(hObject, eventdata, handles)
% hObject    handle to RMA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ѡ��InH�����µĴ����¼�
% ��������ͬ��������ѡ��״̬
set(handles.UMI,'Value',0);
set(handles.UMA,'Value',0);
set(handles.InH,'Value',0);
%����wraparoundѡ�����Ϊ������wraparound
set(handles.wrap,'Value',1);
set(handles.wrap,'Enable','on');
set(handles.Rdrop,'Enable','off');%����UMI����������
set(handles.femto,'Enable','off');%����InH����������
set(handles.slot,'String','5');%���ÿ�����
%���û�վ��ز���
set(handles.p_bs1,'String','43');
set(handles.h_bs1,'String','35');
set(handles.bs_max_gain1,'String','8');
%��ղ�����ƫ��״̬������
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
% % set(handles.loc,'Value',1);
% % set(handles.loc,'Enable','on');
set(handles.res,'String','50'); %����Ƶ�ʵ�����
%����Ƶ�ʵ����á���900��RMA������Ƶ��Ҫ��
set(handles.f1,'String','7');
set(handles.f2,'String','7.2');
%����������ISD����
set(handles.isd,'String','1732');
set(handles.isd,'Enable','on');
set(handles.handover,'Value',1);%����handover״̬Ϊ����3dB handover
if get(handles.sysnum,'Value')==1 %�����˫ϵͳ�������ַ��ѡ������Ĭ��ֵΪ��ַ
    set(handles.loc,'Value',1);
    set(handles.loc,'Enable','on');
end
set(handles.shadow_r,'String','0.5');%������Ӱ˥������ϵ��Ϊ0.5
% ����ϵͳ����ѡ��
set(handles.v_type,'Enable','on');
set(handles.a_type,'Enable','on');

% set(handles.sysnum,'Value',1);
% Hint: get(hObject,'Value') returns toggle state of RMA


% --- Executes on button press in UMA.
function UMA_Callback(hObject, eventdata, handles)
% hObject    handle to UMA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ѡ��InH�����µĴ����¼�
% ��������ͬ��������ѡ��״̬
set(handles.UMI,'Value',0);
set(handles.RMA,'Value',0);
set(handles.InH,'Value',0);
%����wraparoundѡ�����Ϊ������wraparound
set(handles.wrap,'Value',1);
set(handles.wrap,'Enable','on');
set(handles.Rdrop,'Enable','off');%����UMI����������
set(handles.femto,'Enable','off');%����InH����������
set(handles.slot,'String','5');%���ÿ�����
%���û�վ��ز���
set(handles.p_bs1,'String','43');
set(handles.h_bs1,'String','25');
set(handles.bs_max_gain1,'String','8');
%��ղ�����ƫ��״̬�����á������������ù���ת�Ƶ���ַ״̬ѡ��Ĵ����¼���
set(handles.div_angle,'String','');
set(handles.div_dist,'String','');
set(handles.div_angle,'Enable','off');
set(handles.div_dist,'Enable','off');
% set(handles.loc,'Value',1);
% set(handles.loc,'Enable','on');
set(handles.res,'String','5');%���÷ֱ���
%����Ƶ�ʵ�����
set(handles.f1,'String','30');
set(handles.f2,'String','30.2');
%����������ISD����
set(handles.isd,'String','200');
set(handles.isd,'Enable','on');
set(handles.handover,'Value',1); %����handover״̬Ϊ������3dB handover
if get(handles.sysnum,'Value')==1 %�����˫ϵͳ�������ַ��ѡ������Ĭ��ֵΪ��ַ
    set(handles.loc,'Value',1);
    set(handles.loc,'Enable','on');
end
set(handles.shadow_r,'String','0.5'); %������Ӱ˥������ϵ��Ϊ0.5
% ����ϵͳ����ѡ��
set(handles.v_type,'Enable','on');
set(handles.a_type,'Enable','on');
% set(handles.sysnum,'Value',1);

% Hint: get(hObject,'Value') returns toggle state of UMA


% --- Executes on selection change in together.
function together_Callback(hObject, eventdata, handles)
% hObject    handle to together (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns together contents as cell array
%        contents{get(hObject,'Value')} returns selected item from together


% --- Executes during object creation, after setting all properties.
function together_CreateFcn(hObject, eventdata, handles)
% hObject    handle to together (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% ѡ��������ͣ�ͬ��/�첽���Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function div_angle_Callback(hObject, eventdata, handles)
% hObject    handle to div_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ϵͳƫ�ýǶȵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of div_angle as text
%        str2double(get(hObject,'String')) returns contents of div_angle as a double


% --- Executes during object creation, after setting all properties.
function div_angle_CreateFcn(hObject, eventdata, handles)
% hObject    handle to div_angle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function div_dist_Callback(hObject, eventdata, handles)
% hObject    handle to div_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ϵͳƫ�þ���Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of div_dist as text
%        str2double(get(hObject,'String')) returns contents of div_dist as a double


% --- Executes during object creation, after setting all properties.
function div_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to div_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function res_Callback(hObject, eventdata, handles)
% hObject    handle to res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ϵͳ�ֱ��ʵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of res as text
%        str2double(get(hObject,'String')) returns contents of res as a double


% --- Executes during object creation, after setting all properties.
function res_CreateFcn(hObject, eventdata, handles)
% hObject    handle to res (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function isd_Callback(hObject, eventdata, handles)
% hObject    handle to isd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ϵͳվ���Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of isd as text
%        str2double(get(hObject,'String')) returns contents of isd as a double


% --- Executes during object creation, after setting all properties.
function isd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to isd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in wrap.
function wrap_Callback(hObject, eventdata, handles)
% hObject    handle to wrap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ϵͳ�Ƿ����wraparound�Ĵ����¼�
% ��Ҫ�޸ķֱ���
if get(handles.UMA,'Value') && get(handles.wrap,'Value')==2
    set(handles.res,'String','10');
elseif get(handles.UMA,'Value') && get(handles.wrap,'Value')==1
    set(handles.res,'String','5');
elseif get(handles.RMA,'Value')
    set(handles.res,'String','80');
end
    
    
% Hints: contents = cellstr(get(hObject,'String')) returns wrap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from wrap


% --- Executes during object creation, after setting all properties.
function wrap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wrap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f1_Callback(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����victimϵͳƵ�ʵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of f1 as text
%        str2double(get(hObject,'String')) returns contents of f1 as a double


% --- Executes during object creation, after setting all properties.
function f1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function band_Callback(hObject, eventdata, handles)
% hObject    handle to band (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ϵͳ����Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of band as text
%        str2double(get(hObject,'String')) returns contents of band as a double


% --- Executes during object creation, after setting all properties.
function band_CreateFcn(hObject, eventdata, handles)
% hObject    handle to band (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p_ue_Callback(hObject, eventdata, handles)
% hObject    handle to p_ue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ϵͳUE���书�ʵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of p_ue as text
%        str2double(get(hObject,'String')) returns contents of p_ue as a double


% --- Executes during object creation, after setting all properties.
function p_ue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_ue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nf_ue_Callback(hObject, eventdata, handles)
% hObject    handle to nf_ue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ϵͳUE������ϵ���Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of nf_ue as text
%        str2double(get(hObject,'String')) returns contents of nf_ue as a double


% --- Executes during object creation, after setting all properties.
function nf_ue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nf_ue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function v_type_Callback(hObject, eventdata, handles)
% hObject    handle to v_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ��victimϵͳ�����ʹ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of v_type as text
%        str2double(get(hObject,'String')) returns contents of v_type as a double


% --- Executes during object creation, after setting all properties.
function v_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function h_bs1_Callback(hObject, eventdata, handles)
% hObject    handle to h_bs1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����victimϵͳ��ͬ������ϵͳ�Ļ�վ�߶ȵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of h_bs1 as text
%        str2double(get(hObject,'String')) returns contents of h_bs1 as a double


% --- Executes during object creation, after setting all properties.
function h_bs1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h_bs1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function f2_Callback(hObject, eventdata, handles)
% hObject    handle to f2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����aggressorϵͳƵ�ʵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of f2 as text
%        str2double(get(hObject,'String')) returns contents of f2 as a double


% --- Executes during object creation, after setting all properties.
function f2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to f2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p_bs1_Callback(hObject, eventdata, handles)
% hObject    handle to p_bs1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����victimϵͳ��ͬ������ϵͳ�Ļ�վ���书�ʵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of p_bs1 as text
%        str2double(get(hObject,'String')) returns contents of p_bs1 as a double


% --- Executes during object creation, after setting all properties.
function p_bs1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_bs1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nf_bs_Callback(hObject, eventdata, handles)
% hObject    handle to nf_bs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ���û�վ������ϵ���Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of nf_bs as text
%        str2double(get(hObject,'String')) returns contents of nf_bs as a double


% --- Executes during object creation, after setting all properties.
function nf_bs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nf_bs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a_type_Callback(hObject, eventdata, handles)
% hObject    handle to a_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����aggressorϵͳ���͵Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of a_type as text
%        str2double(get(hObject,'String')) returns contents of a_type as a double


% --- Executes during object creation, after setting all properties.
function a_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function shadow_r_Callback(hObject, eventdata, handles)
% hObject    handle to shadow_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ϵͳ��Ӱ˥��ϵ���Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of shadow_r as text
%        str2double(get(hObject,'String')) returns contents of shadow_r as a double


% --- Executes during object creation, after setting all properties.
function shadow_r_CreateFcn(hObject, eventdata, handles)
% hObject    handle to shadow_r (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in handover.
function handover_Callback(hObject, eventdata, handles)
% hObject    handle to handover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ���Ƿ����3dB handover�Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: contents = cellstr(get(hObject,'String')) returns handover contents as cell array
%        contents{get(hObject,'Value')} returns selected item from handover


% --- Executes during object creation, after setting all properties.
function handover_CreateFcn(hObject, eventdata, handles)
% hObject    handle to handover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function p_bs2_Callback(hObject, eventdata, handles)
% hObject    handle to p_bs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% �����칹����ϵͳ�Ļ�վ���书�ʵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of p_bs2 as text
%        str2double(get(hObject,'String')) returns contents of p_bs2 as a double


% --- Executes during object creation, after setting all properties.
function p_bs2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to p_bs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function h_bs2_Callback(hObject, eventdata, handles)
% hObject    handle to h_bs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����aggressorϵͳ��վ�߶ȵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of h_bs2 as text
%        str2double(get(hObject,'String')) returns contents of h_bs2 as a double


% --- Executes during object creation, after setting all properties.
function h_bs2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to h_bs2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in isphoto.
function isphoto_Callback(hObject, eventdata, handles)
% hObject    handle to isphoto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ���Ƿ�����м���ͼ��Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hint: get(hObject,'Value') returns toggle state of isphoto


% --- Executes on button press in issave.
function issave_Callback(hObject, eventdata, handles)
% hObject    handle to issave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ���Ƿ񱣴沨�����;���Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hint: get(hObject,'Value') returns toggle state of issave


% --- Executes on button press in use_cache.
function use_cache_Callback(hObject, eventdata, handles)
% hObject    handle to use_cache (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% �����Ƿ�ʹ�ò������ͻ���Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hint: get(hObject,'Value') returns toggle state of use_cache


% --- Executes on button press in ue_bf.
function ue_bf_Callback(hObject, eventdata, handles)
% hObject    handle to ue_bf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ��UE�Ƿ����BF�����¼�
if get(hObject,'Value') == 1 %������ֻ��UE�в������ͣ���BS�޲������͵����������UEһ��ѡ�������ͣ�BSҲһ��ѡ��BF
    set(handles.bs_bf,'Value',1)
end
% Hint: get(hObject,'Value') returns toggle state of ue_bf


% --- Executes on button press in bs_bf.
function bs_bf_Callback(hObject, eventdata, handles)
% hObject    handle to bs_bf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ��BS�Ƿ����BF�����¼�
if get(hObject,'Value') == 0 %ͬ�ϣ�BSһ��û�в������ͣ�UEһ��û�в�������
    set(handles.ue_bf,'Value',0)
end
% Hint: get(hObject,'Value') returns toggle state of bs_bf



function bs_array_out_Callback(hObject, eventdata, handles)
% hObject    handle to bs_array_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% �������ⳡ��BS������Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of bs_array_out as text
%        str2double(get(hObject,'String')) returns contents of bs_array_out as a double


% --- Executes during object creation, after setting all properties.
function bs_array_out_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_array_out (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ue_max_gain_Callback(hObject, eventdata, handles)
% hObject    handle to ue_max_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����UE�����Ԫ����Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of ue_max_gain as text
%        str2double(get(hObject,'String')) returns contents of ue_max_gain as a double


% --- Executes during object creation, after setting all properties.
function ue_max_gain_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ue_max_gain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bs_phi_Callback(hObject, eventdata, handles)
% hObject    handle to bs_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����BS����ˮƽ���������Ͻ�Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of bs_phi as text
%        str2double(get(hObject,'String')) returns contents of bs_phi as a double


% --- Executes during object creation, after setting all properties.
function bs_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ue_phi_Callback(hObject, eventdata, handles)
% hObject    handle to ue_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����UE����ˮƽ���������Ͻ�Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of ue_phi as text
%        str2double(get(hObject,'String')) returns contents of ue_phi as a double


% --- Executes during object creation, after setting all properties.
function ue_phi_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ue_phi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bs_max_gain2_Callback(hObject, eventdata, handles)
% hObject    handle to bs_max_gain2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% �����칹ϵͳ��վ����书�ʵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of bs_max_gain2 as text
%        str2double(get(hObject,'String')) returns contents of bs_max_gain2 as a double


% --- Executes during object creation, after setting all properties.
function bs_max_gain2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_max_gain2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bs_max_gain1_Callback(hObject, eventdata, handles)
% hObject    handle to bs_max_gain1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����victimϵͳ��ͬ������ϵͳ��վ���书�ʵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of bs_max_gain1 as text
%        str2double(get(hObject,'String')) returns contents of bs_max_gain1 as a double


% --- Executes during object creation, after setting all properties.
function bs_max_gain1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_max_gain1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bs_theta_Callback(hObject, eventdata, handles)
% hObject    handle to bs_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����BS������ֱ������Ͻ�Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of bs_theta as text
%        str2double(get(hObject,'String')) returns contents of bs_theta as a double


% --- Executes during object creation, after setting all properties.
function bs_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ue_theta_Callback(hObject, eventdata, handles)
% hObject    handle to ue_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����UE������ֱ������Ͻ�Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of ue_theta as text
%        str2double(get(hObject,'String')) returns contents of ue_theta as a double


% --- Executes during object creation, after setting all properties.
function ue_theta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ue_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function acir_lower_Callback(hObject, eventdata, handles)
% hObject    handle to acir_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% �����о�ACIR�½�Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of acir_lower as text
%        str2double(get(hObject,'String')) returns contents of acir_lower as a double


% --- Executes during object creation, after setting all properties.
function acir_lower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acir_lower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function acir_step_Callback(hObject, eventdata, handles)
% hObject    handle to acir_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ACIR�о������Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of acir_step as text
%        str2double(get(hObject,'String')) returns contents of acir_step as a double


% --- Executes during object creation, after setting all properties.
function acir_step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acir_step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function slot_Callback(hObject, eventdata, handles)
% hObject    handle to slot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ���ÿ������Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of slot as text
%        str2double(get(hObject,'String')) returns contents of slot as a double


% --- Executes during object creation, after setting all properties.
function slot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function acir_upper_Callback(hObject, eventdata, handles)
% hObject    handle to acir_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% �������о���ACIR�Ͻ�Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of acir_upper as text
%        str2double(get(hObject,'String')) returns contents of acir_upper as a double


% --- Executes during object creation, after setting all properties.
function acir_upper_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acir_upper (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function acir_Callback(hObject, eventdata, handles)
% hObject    handle to acir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ACIRУ׼ֵ�Ĵ����¼�
% ֻ�贫���趨ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of acir as text
%        str2double(get(hObject,'String')) returns contents of acir as a double


% --- Executes during object creation, after setting all properties.
function acir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to acir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in distribution.
function distribution_Callback(hObject, eventdata, handles)
% hObject    handle to distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ�������ֲ��Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: contents = cellstr(get(hObject,'String')) returns distribution contents as cell array
%        contents{get(hObject,'Value')} returns selected item from distribution


% --- Executes during object creation, after setting all properties.
function distribution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to distribution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in still.
function still_Callback(hObject, eventdata, handles)
% hObject    handle to still (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ��UE��ֹ�Ĵ����¼�

set(handles.move,'Value',0); % ��Ҫ����UE�ƶ���ѡ��
set(handles.speed,'Enable','off'); %�����ٶ�����ѡ��
set(handles.TTI,'Enable','off'); %��������ʱ������
% Hint: get(hObject,'Value') returns toggle state of still



function speed_Callback(hObject, eventdata, handles)
% hObject    handle to speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% �����ƶ��ٶȵĴ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of speed as text
%        str2double(get(hObject,'String')) returns contents of speed as a double


% --- Executes during object creation, after setting all properties.
function speed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to speed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TTI_Callback(hObject, eventdata, handles)
% hObject    handle to TTI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ���÷���ʱ���Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of TTI as text
%        str2double(get(hObject,'String')) returns contents of TTI as a double


% --- Executes during object creation, after setting all properties.
function TTI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TTI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in move.
function move_Callback(hObject, eventdata, handles)
% hObject    handle to move (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ѡ��UE�ƶ��Ĵ����¼�
set(handles.still,'Value',0); %����UE��ֹ��ѡ��
set(handles.speed,'Enable','on'); %����UE�ٶ�ѡ��
set(handles.TTI,'Enable','on'); %��������ʱ��
% Hint: get(hObject,'Value') returns toggle state of move



function bs_array_in_Callback(hObject, eventdata, handles)
% hObject    handle to bs_array_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% �������ڻ�վ�������ģ�Ĵ����¼�
% ֻ�贫��ѡ��ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of bs_array_in as text
%        str2double(get(hObject,'String')) returns contents of bs_array_in as a double


% --- Executes during object creation, after setting all properties.
function bs_array_in_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bs_array_in (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function L_2_F_Callback(hObject, eventdata, handles)
% hObject    handle to isphoto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% �칹�����к��߸���ǰ��
set(handles.F_2_L,'Value',0); %ǰ�߸��ź�������߸���ǰ��ֻ��ѡ������һ��
if get(handles.UMA2UMI,'Value')+get(handles.UMA2InH,'Value')>0 %����칹ϵͳ����UMA
    set(handles.v_type,'Enable','on');
    set(handles.a_type,'Enable','off');
    set(handles.v_type,'Value',1);
    set(handles.a_type,'Value',1);
end
% Hint: get(hObject,'Value') returns toggle state of L_2_F

function F_2_L_Callback(hObject, eventdata, handles)
% hObject    handle to isphoto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%�칹������ǰ�߸��ź���
set(handles.L_2_F,'Value',0); %ǰ�߸��ź�������߸���ǰ��ֻ��ѡ������һ��
if get(handles.UMA2UMI,'Value')+get(handles.UMA2InH,'Value')>0 %����칹ϵͳ����UMA
    set(handles.v_type,'Enable','off');
    set(handles.a_type,'Enable','on');
    set(handles.v_type,'Value',1);
    set(handles.a_type,'Value',1);
end
% Hint: get(hObject,'Value') returns toggle state of F_2_L

function UEs_per_slot_Callback(hObject, eventdata, handles)
% hObject    handle to UEs_per_slot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% ����ÿ�ο��շ����UE��
% ֻ�贫���趨ֵ���ɣ�������Լ����ϵ
% Hints: get(hObject,'String') returns contents of UEs_per_slot as text
%        str2double(get(hObject,'String')) returns contents of UEs_per_slot as a double


% --- Executes during object creation, after setting all properties.
function UEs_per_slot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UEs_per_slot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
