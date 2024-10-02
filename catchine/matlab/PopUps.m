function PopUps(s,var)

% pop up menus from structure

names = fieldnames(s);
n = length(names);

posx = 0;
szx = 120;
dx = 10;
posy = 0;
szy = 25;
dy = 10;

for i=1:n
    name = names{i};
    uicontrol('Style','edit',...
        'Position',[posx i*(szy+dy)+7 szx szy-10],...
        'String',[name,'  '],...
        'HorizontalAlignment','right',...
        'Enable','inactive',...
        'Style','text');
    val = getfield(s,name);
    if iscell(val) && ischar(val{1})
        my_callback = @(hObj,event) evalin('caller',[var,'.',name,'=''',val{get(hObj,'Value')},''';']);
        evalin('caller',[var,'.',name,'=''',val{1},''';']);
        uicontrol('Style', 'popup',...
            'String', cell2tag(val),...
            'Position', [posx+szx+dx i*(szy+dy) szx szy],...
            'Callback', my_callback);
    elseif iscell(val) && isnumeric(val{1})
        my_callback = @(hObj,event) evalin('caller',[var,'.',name,'=',num2str(val{get(hObj,'Value')}),';']);
        evalin('caller',[var,'.',name,'=',num2str(val{1}),';']);
        uicontrol('Style', 'popup',...
            'String', cell2tag(val),...
            'Position', [posx+szx+dx i*(szy+dy) szx szy],...
            'Callback', my_callback);
    else
        my_callback = @(hObj,event) evalin('caller',[var,'.',name,'=',get(hObj,'String'),';']);
        evalin('caller',[var,'.',name,'=',num2str(val(1)),';']);
        h=uicontrol('Style', 'edit','String',num2str(val(1)),...
            'Position', [posx+szx+dx i*(szy+dy) szx szy],...
            'HorizontalAlignment','center',...
            'Callback', my_callback);
    end
    
    
end





function tag = cell2tag(c)
tag = [];
for i=1:length(c)
    if ischar(c{i})
        tag = [tag,'|',c{i}];
    else % assume numeric
        tag = [tag,'|',num2str(c{i})];
    end
end
tag = tag(2:end);