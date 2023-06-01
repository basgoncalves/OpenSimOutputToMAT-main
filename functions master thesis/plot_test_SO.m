function plot_test_SO
clc;clear all;

op = '0001';

hip_op_f = '100';
knee_op_f = '100';
% 
% % modify actuator file
% act = xml_read('C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\TD06\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\3DGAIT_B_W14\Output\SO\FINAL_PERSONALISEDTORSIONS_scaled_final_actuators.xml');
% for i = 1:17
%     if contains(act.ForceSet.objects.CoordinateActuator(i).coordinate,'hip_flexion_l')
%         act.ForceSet.objects.CoordinateActuator(i).optimal_force = hip_op_f;
%     elseif contains(act.ForceSet.objects.CoordinateActuator(i).coordinate,'knee_angle_l')
%         act.ForceSet.objects.CoordinateActuator(i).optimal_force = knee_op_f;
%     end
% end
% % save new actuator file
% root = 'OpenSimDocument';                                                        
% Pref = struct;
% Pref.StructItem = false;
% Pref.CellItem = false;
% act = ConvertLogicToString(act);
% path_file_output = ['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\TD06\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\3DGAIT_B_W14\Output\SO\FINAL_PERSONALISEDTORSIONS_scaled_final_actuators_hip_' hip_op_f '_knee_' knee_op_f '.xml'];
% xml_write(path_file_output, act, root,Pref);
% 
% 
% % modify settings file
% set = xml_read('C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\TD06\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\3DGAIT_B_W14\Output\SO\soSettings.xml');

% import data
% forces = load_sto_file(['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\TD06\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\3DGAIT_B_W14\Output\SO\_op' op '_hip' hip_op_f '_knee' knee_op_f '_StaticOptimization_force.sto']);
leg = 'l';
forces = load_sto_file(['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\TD06\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\3DGAIT_B_W14\Output\SO\_StaticOptimization_force.sto']);
t = '14';
% % forces = load_sto_file(['C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\TD06\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\3DGAIT_B_W' t '\Output\SO\_new__StaticOptimization_force.sto']);
switch t
    case '12'
        leg = 'r';
    case '14'
        leg = 'l';
    case '19'
%         start = 7;
%         stop = 103;
        leg = 'r';
    case '20'
%         start = 1;
%         stop = 97;
        leg = 'l';
end

% forces = load_sto_file('C:\Users\Balu\Desktop\_StaticOptimization_force.sto');
med_gas = forces.(['med_gas_' leg])%(start:stop);
psoas = forces.(['psoas_' leg])%(start:stop);
vast_lat = forces.(['vas_lat_' leg])%(start:stop);
rect_fem = forces.(['rect_fem_' leg])%(start:stop);
glut_max2 = forces.(['glut_max2_' leg])%(start:stop);

% plot data
figure
% plot(forces.time(start:stop),med_gas)
% hold on
% plot(forces.time(start:stop),psoas)
% plot(forces.time(start:stop),vast_lat)
% plot(forces.time(start:stop),rect_fem)
% plot(forces.time(start:stop),glut_max2)
plot(forces.time,med_gas)
hold on
plot(forces.time,psoas)
plot(forces.time,vast_lat)
plot(forces.time,rect_fem)
plot(forces.time,glut_max2)
% title(['muscle forces: hip ' hip_op_f ', knee ' knee_op_f ' and op 0.' op 'no settings file'])
title(['muscle forces Trial ' t])

legend('gastrocnemius medialis', 'psoas', 'vastus lateralis', 'rectus femoris', 'gluteus maximus 2')
xlabel('time [s]')
ylabel('muscle force [N]')

%save data
savepath = 'C:\Users\Balu\Nextcloud\Documents\MA\Code\Kira_MSc_data\TD06\output automization\FINAL_PERSONALISEDTORSIONS_scaled_final\3DGAIT_B_W14\Output\SO';
% saveas(gcf,fullfile(savepath, ['muscle_forces_l_op' op '_hip' hip_op_f '_knee' knee_op_f]),'tif')





% ------------------------------------------------------------------- %
% -------------------------- FUNCTIONS -------------------------------%

function [tree, RootName, DOMnode] = xml_read(xmlfile, Pref)
%XML_READ reads xml files and converts them into Matlab's struct tree.
%
% DESCRIPTION
% tree = xml_read(xmlfile) reads 'xmlfile' into data structure 'tree'
%
% tree = xml_read(xmlfile, Pref) reads 'xmlfile' into data structure 'tree'
% according to your preferences
%
% [tree, RootName, DOMnode] = xml_read(xmlfile) get additional information
% about XML file
%
% INPUT:
%  xmlfile	URL or filename of xml file to read
%  Pref     Preferences:
%    Pref.ItemName - default 'item' - name of a special tag used to itemize
%                    cell arrays
%    Pref.ReadAttr - default true - allow reading attributes
%    Pref.ReadSpec - default true - allow reading special nodes
%    Pref.Str2Num  - default 'smart' - convert strings that look like numbers
%                   to numbers. Options: "always", "never", and "smart"
%    Pref.KeepNS   - default true - keep or strip namespace info
%    Pref.NoCells  - default true - force output to have no cell arrays
%    Pref.Debug    - default false - show mode specific error messages
%    Pref.NumLevels- default infinity - how many recursive levels are
%      allowed. Can be used to speed up the function by prunning the tree.
%    Pref.RootOnly - default true - output variable 'tree' corresponds to
%      xml file root element, otherwise it correspond to the whole file.
%    Pref.CellItem - default 'true' - leave 'item' nodes in cell notation.
% OUTPUT:
%  tree         tree of structs and/or cell arrays corresponding to xml file
%  RootName     XML tag name used for root (top level) node.
%               Optionally it can be a string cell array storing: Name of
%               root node, document "Processing Instructions" data and
%               document "comment" string
%  DOMnode      output of xmlread
%
% DETAILS:
% Function xml_read first calls MATLAB's xmlread function and than
% converts its output ('Document Object Model' tree of Java objects)
% to tree of MATLAB struct's. The output is in format of nested structs
% and cells. In the output data structure field names are based on
% XML tags, except in cases when tags produce illegal variable names.
%
% Several special xml node types result in special tags for fields of
% 'tree' nodes:
%  - node.CONTENT - stores data section of the node if other fields are
%    present. Usually data section is stored directly in 'node'.
%  - node.ATTRIBUTE.name - stores node's attribute called 'name'.
%  - node.COMMENT - stores node's comment section (string). For global
%    comments see "RootName" output variable.
%  - node.CDATA_SECTION - stores node's CDATA section (string).
%  - node.PROCESSING_INSTRUCTIONS - stores "processing instruction" child
%    node. For global "processing instructions" see "RootName" output variable.
%  - other special node types like: document fragment nodes, document type
%   nodes, entity nodes, notation nodes and processing instruction nodes
%   will be treated like regular nodes
%
% EXAMPLES:
%   MyTree=[];
%   MyTree.MyNumber = 13;
%   MyTree.MyString = 'Hello World';
%   xml_write('test.xml', MyTree);
%   [tree treeName] = xml_read ('test.xml');
%   disp(treeName)
%   gen_object_display()
%   % See also xml_examples.m
%
% See also:
%   xml_write, xmlread, xmlwrite
%
% Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com
% References:
%  - Function inspired by Example 3 found in xmlread function.
%  - Output data structures inspired by xml_toolbox structures.

%% default preferences
DPref.TableName  = {'tr','td'}; % name of a special tags used to itemize 2D cell arrays
DPref.ItemName  = 'item'; % name of a special tag used to itemize 1D cell arrays
DPref.CellItem  = false;  % leave 'item' nodes in cell notation
DPref.ReadAttr  = true;   % allow reading attributes
DPref.ReadSpec  = true;   % allow reading special nodes: comments, CData, etc.
DPref.KeepNS    = true;   % Keep or strip namespace info
DPref.Str2Num   = 'smart';% convert strings that look like numbers to numbers
DPref.NoCells   = true;   % force output to have no cell arrays
DPref.NumLevels = 1e10;   % number of recurence levels
DPref.PreserveSpace = false; % Preserve or delete spaces at the beggining and the end of stings?
RootOnly        = true;   % return root node  with no top level special nodes
Debug           = false;  % show specific errors (true) or general (false)?
tree            = [];
RootName        = [];

%% Check Matlab Version
v = ver('MATLAB');
version = str2double(regexp(v.Version, '\d.\d','match','once'));
if (version<7.1)
  error('Your MATLAB version is too old. You need version 7.1 or newer.');
end

%% read user preferences
if (nargin>1)
  if (isfield(Pref, 'TableName')), DPref.TableName = Pref.TableName; end
  if (isfield(Pref, 'ItemName' )), DPref.ItemName  = Pref.ItemName;  end
  if (isfield(Pref, 'CellItem' )), DPref.CellItem  = Pref.CellItem;  end
  if (isfield(Pref, 'Str2Num'  )), DPref.Str2Num   = Pref.Str2Num ;  end
  if (isfield(Pref, 'NoCells'  )), DPref.NoCells   = Pref.NoCells ;  end
  if (isfield(Pref, 'NumLevels')), DPref.NumLevels = Pref.NumLevels; end
  if (isfield(Pref, 'ReadAttr' )), DPref.ReadAttr  = Pref.ReadAttr;  end
  if (isfield(Pref, 'ReadSpec' )), DPref.ReadSpec  = Pref.ReadSpec;  end
  if (isfield(Pref, 'KeepNS'   )), DPref.KeepNS    = Pref.KeepNS;    end
  if (isfield(Pref, 'RootOnly' )), RootOnly        = Pref.RootOnly;  end
  if (isfield(Pref, 'Debug'    )), Debug           = Pref.Debug   ;  end
  if (isfield(Pref, 'PreserveSpace')), DPref.PreserveSpace = Pref.PreserveSpace; end
end
if ischar(DPref.Str2Num), % convert from character description to numbers
  DPref.Str2Num = find(strcmpi(DPref.Str2Num, {'never', 'smart', 'always'}))-1;
  if isempty(DPref.Str2Num), DPref.Str2Num=1; end % 1-smart by default
end

%% read xml file using Matlab function
if isa(xmlfile, 'org.apache.xerces.dom.DeferredDocumentImpl');
  % if xmlfile is a DOMnode than skip the call to xmlread
  try
    try
      DOMnode = xmlfile;
    catch ME
      error('Invalid DOM node: \n%s.', getReport(ME));
    end
  catch %#ok<CTCH> catch for mablab versions prior to 7.5
    error('Invalid DOM node. \n');
  end
else         % we assume xmlfile is a filename
  if (Debug) % in debuging mode crashes are allowed
    DOMnode = xmlread(xmlfile);
  else       % in normal mode crashes are not allowed
    try
      try
        DOMnode = xmlread(xmlfile);
      catch ME
        error('Failed to read XML file %s: \n%s',xmlfile, getReport(ME));
      end
    catch %#ok<CTCH> catch for mablab versions prior to 7.5
      error('Failed to read XML file %s\n',xmlfile);
    end
  end
end
Node = DOMnode.getFirstChild;

%% Find the Root node. Also store data from Global Comment and Processing
%  Instruction nodes, if any.
GlobalTextNodes = cell(1,3);
GlobalProcInst  = [];
GlobalComment   = [];
GlobalDocType   = [];
while (~isempty(Node))
  if (Node.getNodeType==Node.ELEMENT_NODE)
    RootNode=Node;
  elseif (Node.getNodeType==Node.PROCESSING_INSTRUCTION_NODE)
    data   = strtrim(char(Node.getData));
    target = strtrim(char(Node.getTarget));
    GlobalProcInst = [target, ' ', data];
    GlobalTextNodes{2} = GlobalProcInst;
  elseif (Node.getNodeType==Node.COMMENT_NODE)
    GlobalComment = strtrim(char(Node.getData));
    GlobalTextNodes{3} = GlobalComment;
    %   elseif (Node.getNodeType==Node.DOCUMENT_TYPE_NODE)
    %     GlobalTextNodes{4} = GlobalDocType;
  end
  Node = Node.getNextSibling;
end

%% parse xml file through calls to recursive DOMnode2struct function
if (Debug)   % in debuging mode crashes are allowed
  [tree RootName] = DOMnode2struct(RootNode, DPref, 1);
else         % in normal mode crashes are not allowed
  try
    try
      [tree RootName] = DOMnode2struct(RootNode, DPref, 1);
    catch ME
      error('Unable to parse XML file %s: \n %s.',xmlfile, getReport(ME));
    end
  catch %#ok<CTCH> catch for mablab versions prior to 7.5
    error('Unable to parse XML file %s.',xmlfile);
  end
end

%% If there were any Global Text nodes than return them
if (~RootOnly)
  if (~isempty(GlobalProcInst) && DPref.ReadSpec)
    t.PROCESSING_INSTRUCTION = GlobalProcInst;
  end
  if (~isempty(GlobalComment) && DPref.ReadSpec)
    t.COMMENT = GlobalComment;
  end
  if (~isempty(GlobalDocType) && DPref.ReadSpec)
    t.DOCUMENT_TYPE = GlobalDocType;
  end
  t.(RootName) = tree;
  tree=t;
end
if (~isempty(GlobalTextNodes))
  GlobalTextNodes{1} = RootName;
  RootName = GlobalTextNodes;
end


%% =======================================================================
%  === DOMnode2struct Function ===========================================
%  =======================================================================
function [s TagName LeafNode] = DOMnode2struct(node, Pref, level)

%% === Step 1: Get node name and check if it is a leaf node ==============
[TagName LeafNode] = NodeName(node, Pref.KeepNS);
s = []; % initialize output structure

%% === Step 2: Process Leaf Nodes (nodes with no children) ===============
if (LeafNode)
  if (LeafNode>1 && ~Pref.ReadSpec), LeafNode=-1; end % tags only so ignore special nodes
  if (LeafNode>0) % supported leaf node types
    try
      try         % use try-catch: errors here are often due to VERY large fields (like images) that overflow java memory
        s = char(node.getData);
        if (isempty(s)), s = ' '; end                              % make it a string
        % for some reason current xmlread 'creates' a lot of empty text
        % fields with first chatacter=10 - those will be deleted.
        if (~Pref.PreserveSpace || s(1)==10) 
          if (isspace(s(1)) || isspace(s(end))), s = strtrim(s); end % trim speces is any
        end
        if (LeafNode==1), s=str2var(s, Pref.Str2Num, 0); end       % convert to number(s) if needed
      catch ME    % catch for mablab versions 7.5 and higher
        warning('xml_io_tools:read:LeafRead', ...
          'This leaf node could not be read and was ignored. ');
        getReport(ME)
      end
    catch         %#ok<CTCH> catch for mablab versions prior to 7.5
      warning('xml_io_tools:read:LeafRead', ...
        'This leaf node could not be read and was ignored. ');
    end
  end
  if (LeafNode==3) % ProcessingInstructions need special treatment
    target = strtrim(char(node.getTarget));
    s = [target, ' ', s];
  end
  return % We are done the rest of the function deals with nodes with children
end
if (level>Pref.NumLevels+1), return; end % if Pref.NumLevels is reached than we are done

%% === Step 3: Process nodes with children ===============================
if (node.hasChildNodes)        % children present
  Child  = node.getChildNodes; % create array of children nodes
  nChild = Child.getLength;    % number of children
  
  % --- pass 1: how many children with each name -----------------------
  f = [];
  for iChild = 1:nChild        % read in each child
    [cname cLeaf] = NodeName(Child.item(iChild-1), Pref.KeepNS);
    if (cLeaf<0), continue; end % unsupported leaf node types
    if (~isfield(f,cname)),
      f.(cname)=0;           % initialize first time I see this name
    end
    f.(cname) = f.(cname)+1; % add to the counter
  end                        % end for iChild
  % text_nodes become CONTENT & for some reason current xmlread 'creates' a
  % lot of empty text fields so f.CONTENT value should not be trusted
  if (isfield(f,'CONTENT') && f.CONTENT>2), f.CONTENT=2; end
  
  % --- pass 2: store all the children as struct of cell arrays ----------
  for iChild = 1:nChild        % read in each child
    [c cname cLeaf] = DOMnode2struct(Child.item(iChild-1), Pref, level+1);
    if (cLeaf && isempty(c))   % if empty leaf node than skip
      continue;                % usually empty text node or one of unhandled node types
    elseif (nChild==1 && cLeaf==1)
      s=c;                     % shortcut for a common case
    else                       % if normal node
      if (level>Pref.NumLevels), continue; end
      n = f.(cname);           % how many of them in the array so far?
      if (~isfield(s,cname))   % encountered this name for the first time
        if (n==1)              % if there will be only one of them ...
          s.(cname) = c;       % than save it in format it came in
        else                   % if there will be many of them ...
          s.(cname) = cell(1,n);
          s.(cname){1} = c;    % than save as cell array
        end
        f.(cname) = 1;         % initialize the counter
      else                     % already have seen this name
        s.(cname){n+1} = c;    % add to the array
        f.(cname) = n+1;       % add to the array counter
      end
    end
  end   % for iChild
end % end if (node.hasChildNodes)

%% === Step 4: Post-process struct's created for nodes with children =====
if (isstruct(s))
  fields = fieldnames(s);
  nField = length(fields);
  
  % Detect structure that looks like Html table and store it in cell Matrix
  if (nField==1 && strcmpi(fields{1},Pref.TableName{1}))
    tr = s.(Pref.TableName{1});
    fields2 = fieldnames(tr{1});
    if (length(fields2)==1 && strcmpi(fields2{1},Pref.TableName{2}))
      % This seems to be a special structure such that for 
      % Pref.TableName = {'tr','td'} 's' corresponds to 
      %    <tr> <td>M11</td> <td>M12</td> </tr>
      %    <tr> <td>M12</td> <td>M22</td> </tr>
      % Recognize it as encoding for 2D struct
      nr = length(tr);
      for r = 1:nr
        row = tr{r}.(Pref.TableName{2});
        Table(r,1:length(row)) = row; %#ok<AGROW>
      end
      s = Table;
    end
  end

  % --- Post-processing: convert 'struct of cell-arrays' to 'array of structs'
  % Example: let say s has 3 fields s.a, s.b & s.c  and each field is an
  % cell-array with more than one cell-element and all 3 have the same length.
  % Then change it to array of structs, each with single cell.
  % This way element s.a{1} will be now accessed through s(1).a
  vec = zeros(size(fields));
  for i=1:nField, vec(i) = f.(fields{i}); end
  if (numel(vec)>1 && vec(1)>1 && var(vec)==0)  % convert from struct of
    s = cell2struct(struct2cell(s), fields, 1); % arrays to array of struct
  end % if anyone knows better way to do above conversion please let me know.

end

%% === Step 5: Process nodes with attributes =============================
if (node.hasAttributes && Pref.ReadAttr)
  if (~isstruct(s)),              % make into struct if is not already
    ss.CONTENT=s;
    s=ss;
  end
  Attr  = node.getAttributes;     % list of all attributes
  for iAttr = 1:Attr.getLength    % for each attribute
    name  = char(Attr.item(iAttr-1).getName);  % attribute name
    name  = str2varName(name, Pref.KeepNS);    % fix name if needed
    value = char(Attr.item(iAttr-1).getValue); % attribute value
    value = str2var(value, Pref.Str2Num, 1);   % convert to number if possible
    s.ATTRIBUTE.(name) = value;   % save again
  end                             % end iAttr loop
end % done with attributes
if (~isstruct(s)), return; end %The rest of the code deals with struct's

%% === Post-processing: fields of "s"
% convert  'cell-array of structs' to 'arrays of structs'
fields = fieldnames(s);     % get field names
nField = length(fields);
for iItem=1:length(s)       % for each struct in the array - usually one
  for iField=1:length(fields)
    field = fields{iField}; % get field name
    % if this is an 'item' field and user want to leave those as cells
    % than skip this one
    if (strcmpi(field, Pref.ItemName) && Pref.CellItem), continue; end
    x = s(iItem).(field);
    if (iscell(x) && all(cellfun(@isstruct,x(:))) && numel(x)>1) % it's cell-array of structs
      % numel(x)>1 check is to keep 1 cell-arrays created when Pref.CellItem=1
      try                           % this operation fails sometimes
        % example: change s(1).a{1}.b='jack'; s(1).a{2}.b='john'; to
        % more convinient s(1).a(1).b='jack'; s(1).a(2).b='john';
        s(iItem).(field) = [x{:}]';  %#ok<AGROW> % converted to arrays of structs
      catch %#ok<CTCH>
        % above operation will fail if s(1).a{1} and s(1).a{2} have
        % different fields. If desired, function forceCell2Struct can force
        % them to the same field structure by adding empty fields.
        if (Pref.NoCells)
          s(iItem).(field) = forceCell2Struct(x); %#ok<AGROW>
        end
      end % end catch
    end
  end
end

%% === Step 4: Post-process struct's created for nodes with children =====

% --- Post-processing: remove special 'item' tags ---------------------
% many xml writes (including xml_write) use a special keyword to mark
% arrays of nodes (see xml_write for examples). The code below converts
% s.item to s.CONTENT
ItemContent = false;
if (isfield(s,Pref.ItemName))
  s.CONTENT = s.(Pref.ItemName);
  s = rmfield(s,Pref.ItemName);
  ItemContent = Pref.CellItem; % if CellItem than keep s.CONTENT as cells
end

% --- Post-processing: clean up CONTENT tags ---------------------
% if s.CONTENT is a cell-array with empty elements at the end than trim
% the length of this cell-array. Also if s.CONTENT is the only field than
% remove .CONTENT part and store it as s.
if (isfield(s,'CONTENT'))
  if (iscell(s.CONTENT) && isvector(s.CONTENT))
    x = s.CONTENT;
    for i=numel(x):-1:1, if ~isempty(x{i}), break; end; end
    if (i==1 && ~ItemContent)
      s.CONTENT = x{1};   % delete cell structure
    else
      s.CONTENT = x(1:i); % delete empty cells
    end
  end
  if (nField==1)
    if (ItemContent)
      ss = s.CONTENT;       % only child: remove a level but ensure output is a cell-array
      s=[]; s{1}=ss;
    else
      s = s.CONTENT;        % only child: remove a level
    end
  end
end



%% =======================================================================
%  === forceCell2Struct Function =========================================
%  =======================================================================
function s = forceCell2Struct(x)
% Convert cell-array of structs, where not all of structs have the same
% fields, to a single array of structs

%% Convert 1D cell array of structs to 2D cell array, where each row
% represents item in original array and each column corresponds to a unique
% field name. Array "AllFields" store fieldnames for each column
AllFields = fieldnames(x{1});     % get field names of the first struct
CellMat = cell(length(x), length(AllFields));
for iItem=1:length(x)
  fields = fieldnames(x{iItem});  % get field names of the next struct
  for iField=1:length(fields)     % inspect all fieldnames and find those
    field = fields{iField};       % get field name
    col = find(strcmp(field,AllFields),1);
    if isempty(col)               % no column for such fieldname yet
      AllFields = [AllFields; field]; %#ok<AGROW>
      col = length(AllFields);    % create a new column for it
    end
    CellMat{iItem,col} = x{iItem}.(field); % store rearanged data
  end
end
%% Convert 2D cell array to array of structs
s = cell2struct(CellMat, AllFields, 2);

%% =======================================================================
%  === str2var Function ==================================================
%  =======================================================================
function val=str2var(str, option, attribute)
% Can this string 'str' be converted to a number? if so than do it.
val = str;
len = numel(str);
if (len==0    || option==0), return; end % Str2Num="never" of empty string -> do not do enything
if (len>10000 && option==1), return; end % Str2Num="smart" and string is very long -> probably base64 encoded binary
digits = '(Inf)|(NaN)|(pi)|[\t\n\d\+\-\*\.ei EI\[\]\;\,]';
s = regexprep(str, digits, ''); % remove all the digits and other allowed characters
if (~all(~isempty(s)))          % if nothing left than this is probably a number
  if (~isempty(strfind(str, ' '))), option=2; end %if str has white-spaces assume by default that it is not a date string
  if (~isempty(strfind(str, '['))), option=2; end % same with brackets
  str(strfind(str, '\n')) = ';';% parse data tables into 2D arrays, if any
  if (option==1)                % the 'smart' option
    try                         % try to convert to a date, like 2007-12-05
      datenum(str);             % if successful than leave it as string
    catch                       %#ok<CTCH> % if this is not a date than ...
      option=2;                 % ... try converting to a number
    end
  end
  if (option==2)
    if (attribute)
      num = str2double(str);      % try converting to a single number using sscanf function
      if isnan(num), return; end  % So, it wasn't really a number after all    
    else
      num = str2num(str);         %#ok<ST2NM> % try converting to a single number or array using eval function
    end
    if(isnumeric(num) && numel(num)>0), val=num; end % if convertion to a single was succesful than save
  end
elseif ((str(1)=='[' && str(end)==']') || (str(1)=='{' && str(end)=='}')) % this looks like a (cell) array encoded as a string
  try 
    val = eval(str); 
  catch              %#ok<CTCH>
    val = str; 
  end                     
elseif (~attribute)   % see if it is a boolean array with no [] brackets
  str1 = lower(str);
  str1 = strrep(str1, 'false', '0');
  str1 = strrep(str1, 'true' , '1');
  s = regexprep(str1, '[01 \;\,]', ''); % remove all 0/1, spaces, commas and semicolons 
  if (~all(~isempty(s)))          % if nothing left than this is probably a boolean array
    num  = str2num(str1); %#ok<ST2NM>
    if(isnumeric(num) && numel(num)>0), val = (num>0);  end % if convertion was succesful than save as logical
  end
end


%% =======================================================================
%  === str2varName Function ==============================================
%  =======================================================================
function str = str2varName(str, KeepNS)
% convert a sting to a valid matlab variable name
if(KeepNS)
  str = regexprep(str,':','_COLON_', 'once', 'ignorecase');
else
  k = strfind(str,':');
  if (~isempty(k))
    str = str(k+1:end);
  end
end
str = regexprep(str,'-','_DASH_'  ,'once', 'ignorecase');
if (~isvarname(str)) && (~iskeyword(str))
  str = genvarname(str);
end

%% =======================================================================
%  === NodeName Function =================================================
%  =======================================================================


function [Name LeafNode] = NodeName(node, KeepNS)
% get node name and make sure it is a valid variable name in Matlab.
% also get node type:
%   LeafNode=0 - normal element node,
%   LeafNode=1 - text node
%   LeafNode=2 - supported non-text leaf node,
%   LeafNode=3 - supported processing instructions leaf node,
%   LeafNode=-1 - unsupported non-text leaf node
switch (node.getNodeType)
  case node.ELEMENT_NODE
    Name = char(node.getNodeName);% capture name of the node
    Name = str2varName(Name, KeepNS);     % if Name is not a good variable name - fix it
    LeafNode = 0;
  case node.TEXT_NODE
    Name = 'CONTENT';
    LeafNode = 1;
  case node.COMMENT_NODE
    Name = 'COMMENT';
    LeafNode = 2;
  case node.CDATA_SECTION_NODE
    Name = 'CDATA_SECTION';
    LeafNode = 2;
  case node.DOCUMENT_TYPE_NODE
    Name = 'DOCUMENT_TYPE';
    LeafNode = 2;
  case node.PROCESSING_INSTRUCTION_NODE
    Name = 'PROCESSING_INSTRUCTION';
    LeafNode = 3;
  otherwise
    NodeType = {'ELEMENT','ATTRIBUTE','TEXT','CDATA_SECTION', ...
      'ENTITY_REFERENCE', 'ENTITY', 'PROCESSING_INSTRUCTION', 'COMMENT',...
      'DOCUMENT', 'DOCUMENT_TYPE', 'DOCUMENT_FRAGMENT', 'NOTATION'};
    Name = char(node.getNodeName);% capture name of the node
    warning('xml_io_tools:read:unkNode', ...
      'Unknown node type encountered: %s_NODE (%s)', NodeType{node.getNodeType}, Name);
    LeafNode = -1;
end

% ------------------------ ConvertLogicToString ----------------------------%

function S = ConvertLogicToString(S,fieldName)
%% This function is need so the xml_write does not save values of 'true' as [true] logicals which will cause problems with running the scale tool in OpenSim.
if nargin < 2
    fieldName = 'apply';
end

if isstruct(S)                                                                                                      % recursive loop that checks all "apply" fields and changed them to a string
    F = fields(S);
    S = editApplyField(S,F,fieldName);
    for i = 1:length(F)
        s = S.(F{i});
        for ii = 1:length(s)
            s(ii) = ConvertLogicToString(s(ii),fieldName);
        end
        S.(F{i}) = s;
    end
end

function s = editApplyField(s,f,fieldName)
%% converts "apply" flieds from [true] / [false] to 'true' / 'false'

try s.(fieldName);
    if islogical(s.(fieldName)) && length(s.(fieldName))<2
        if s.(fieldName)== 1
            s.(fieldName) = 'true';
        elseif s.(fieldName) == 0
            s.(fieldName) = 'false';
        end
    else
        s.(fieldName) = num2str(s.(fieldName));
    end
catch
end


% -------------------- xml_write ----------------------------------------%
function DOMnode = xml_write(filename, tree, RootName, Pref)
%XML_WRITE  Writes Matlab data structures to XML file
%
% DESCRIPTION
% xml_write( filename, tree) Converts Matlab data structure 'tree' containing
% cells, structs, numbers and strings to Document Object Model (DOM) node
% tree, then saves it to XML file 'filename' using Matlab's xmlwrite
% function. Optionally one can also use alternative version of xmlwrite
% function which directly calls JAVA functions for XML writing without
% MATLAB middleware. This function is provided as a patch to existing
% bugs in xmlwrite (in R2006b).
%
% xml_write(filename, tree, RootName, Pref) allows you to specify
% additional preferences about file format
%
% DOMnode = xml_write([], tree) same as above except that DOM node is
% not saved to the file but returned.
%
% INPUT
%   filename     file name
%   tree         Matlab structure tree to store in xml file.
%   RootName     String with XML tag name used for root (top level) node
%                Optionally it can be a string cell array storing: Name of
%                root node, document "Processing Instructions" data and
%                document "comment" stringrotated.vtp
%   Pref         Other preferences:
%     Pref.ItemName - default 'item' -  name of a special tag used to
%                     itemize cell or struct arrays
%     Pref.XmlEngine - let you choose the XML engine. Currently default is
%       'Xerces', which is using directly the apache xerces java file.
%       Other option is 'Matlab' which uses MATLAB's xmlwrite and its
%       XMLUtils java file. Both options create identical results except in
%       case of CDATA sections where xmlwrite fails.
%     Pref.CellItem - default 'true' - allow cell arrays to use 'item'
%       notation. See below.
%    Pref.RootOnly - default true - output variable 'tree' corresponds to
%       xml file root element, otherwise it correspond to the whole file.
%     Pref.StructItem - default 'true' - allow arrays of structs to use
%       'item' notation. For example "Pref.StructItem = true" gives:
%         <a>
%           <b>
%             <item> ... <\item>
%             <item> ... <\item>
%           <\b>
%         <\a>
%       while "Pref.StructItem = false" gives:
%         <a>
%           <b> ... <\b>
%           <b> ... <\b>
%         <\a>
%
%
% Several special xml node types can be created if special tags are used
% for field names of 'tree' nodes:
%  - node.CONTENT - stores data section of the node if other fields
%    (usually ATTRIBUTE are present. Usually data section is stored
%    directly in 'node'.
%  - node.ATTRIBUTE.name - stores node's attribute called 'name'.
%  - node.COMMENT - create comment child node from the string. For global
%    comments see "RootName" input variable.
%  - node.PROCESSING_INSTRUCTIONS - create "processing instruction" child
%    node from the string. For global "processing instructions" see
%    "RootName" input variable.
%  - node.CDATA_SECTION - stores node's CDATA section (string). Only works
%    if Pref.XmlEngine='Xerces'. For more info, see comments of F_xmlwrite.
%  - other special node types like: document fragment nodes, document type
%    nodes, entity nodes and notation nodes are not being handled by
%    'xml_write' at the moment.
%
% OUTPUT
%   DOMnode      Document Object Model (DOM) node tree in the format
%                required as input to xmlwrite. (optional)
%
% EXAMPLES:
%   MyTree=[];
%   MyTree.MyNumber = 13;
%   MyTree.MyString = 'Hello World';
%   xml_write('test.xml', MyTree);
%   type('test.xml')
%   %See also xml_tutorial.m
%
% See also
%   xml_read, xmlread, xmlwrite
%
% Written by Jarek Tuszynski, SAIC, jaroslaw.w.tuszynski_at_saic.com

%% Check Matlab Version
v = ver('MATLAB');
v = str2double(regexp(v.Version, '\d.\d','match','once'));
if (v<7)
    error('Your MATLAB version is too old. You need version 7.0 or newer.');
end

%% default preferences
DPref.TableName  = {'tr','td'}; % name of a special tags used to itemize 2D cell arrays
DPref.ItemName   = 'item'; % name of a special tag used to itemize 1D cell arrays
DPref.StructItem = true;  % allow arrays of structs to use 'item' notation
DPref.CellItem   = true;  % allow cell arrays to use 'item' notation
DPref.StructTable= 'Html';
DPref.CellTable  = 'Html';
DPref.XmlEngine  = 'Matlab';  % use matlab provided XMLUtils
%DPref.XmlEngine  = 'Xerces';  % use Xerces xml generator directly
DPref.PreserveSpace = false; % Preserve or delete spaces at the beggining and the end of stings?
RootOnly         = true;  % Input is root node only
GlobalProcInst = [];
GlobalComment  = [];
GlobalDocType  = [];

%% read user preferences
if (nargin>3)
    if (isfield(Pref, 'TableName' )),  DPref.TableName  = Pref.TableName; end
    if (isfield(Pref, 'ItemName'  )), DPref.ItemName   = Pref.ItemName;   end
    if (isfield(Pref, 'StructItem')), DPref.StructItem = Pref.StructItem; end
    if (isfield(Pref, 'CellItem'  )), DPref.CellItem   = Pref.CellItem;   end
    if (isfield(Pref, 'CellTable')),   DPref.CellTable  = Pref.CellTable; end
    if (isfield(Pref, 'StructTable')), DPref.StructTable= Pref.StructTable; end
    if (isfield(Pref, 'XmlEngine' )), DPref.XmlEngine  = Pref.XmlEngine;  end
    if (isfield(Pref, 'RootOnly'  )), RootOnly         = Pref.RootOnly;   end
    if (isfield(Pref, 'PreserveSpace')), DPref.PreserveSpace = Pref.PreserveSpace; end
else
    Pref.StructItem = false;
end
if (nargin<3 || isempty(RootName)), RootName=inputname(2); end
if (isempty(RootName)), RootName='ROOT'; end
if (iscell(RootName)) % RootName also stores global text node data
    rName = RootName;
    RootName = char(rName{1});
    if (length(rName)>1), GlobalProcInst = char(rName{2}); end
    if (length(rName)>2), GlobalComment  = char(rName{3}); end
    if (length(rName)>3), GlobalDocType  = char(rName{4}); end
end
if(~RootOnly && isstruct(tree))  % if struct than deal with each field separatly
    fields = fieldnames(tree);
    for i=1:length(fields)
        field = fields{i};
        x = tree(1).(field);
        if (strcmp(field, 'COMMENT'))
            GlobalComment = x;
        elseif (strcmp(field, 'PROCESSING_INSTRUCTION'))
            GlobalProcInst = x;
        elseif (strcmp(field, 'DOCUMENT_TYPE'))
            GlobalDocType = x;
        else
            RootName = field;
            t = x;
        end
    end
    tree = t;
end

%% Initialize jave object that will store xml data structure
RootName = varName2str(RootName);
if (~isempty(GlobalDocType))
    %   n = strfind(GlobalDocType, ' ');
    %   if (~isempty(n))
    %     dtype = com.mathworks.xml.XMLUtils.createDocumentType(GlobalDocType);
    %   end
    %   DOMnode = com.mathworks.xml.XMLUtils.createDocument(RootName, dtype);
    warning('xml_io_tools:write:docType', ...
        'DOCUMENT_TYPE node was encountered which is not supported yet. Ignoring.');
end
DOMnode = com.mathworks.xml.XMLUtils.createDocument(RootName);


%% Use recursive function to convert matlab data structure to XML
root = DOMnode.getDocumentElement;
struct2DOMnode(DOMnode, root, tree, DPref.ItemName, DPref);

%% Remove the only child of the root node
root   = DOMnode.getDocumentElement;
Child  = root.getChildNodes; % create array of children nodes
nChild = Child.getLength;    % number of children
if (nChild==1)
    node = root.removeChild(root.getFirstChild);
    while(node.hasChildNodes)
        root.appendChild(node.removeChild(node.getFirstChild));
    end
    while(node.hasAttributes)            % copy all attributes
        root.setAttributeNode(node.removeAttributeNode(node.getAttributes.item(0)));
    end
end

%% Save exotic Global nodes
if (~isempty(GlobalComment))
    DOMnode.insertBefore(DOMnode.createComment(GlobalComment), DOMnode.getFirstChild());
end
if (~isempty(GlobalProcInst))
    n = strfind(GlobalProcInst, ' ');
    if (~isempty(n))
        proc = DOMnode.createProcessingInstruction(GlobalProcInst(1:(n(1)-1)),...
            GlobalProcInst((n(1)+1):end));
        DOMnode.insertBefore(proc, DOMnode.getFirstChild());
    end
end
% Not supported yet as the code below does not work
% if (~isempty(GlobalDocType))
%   n = strfind(GlobalDocType, ' ');
%   if (~isempty(n))
%     dtype = DOMnode.createDocumentType(GlobalDocType);
%     DOMnode.insertBefore(dtype, DOMnode.getFirstChild());
%   end
% end

%% save java DOM tree to XML file
if (~isempty(filename))
    if (strcmpi(DPref.XmlEngine, 'Xerces'))
        xmlwrite_xerces(filename, DOMnode);
    else
        xmlwrite(filename, DOMnode);
    end
end


%% =======================================================================
%  === struct2DOMnode Function ===========================================
%  =======================================================================
function [] = struct2DOMnode(xml, parent, s, TagName, Pref)
% struct2DOMnode is a recursive function that converts matlab's structs to
% DOM nodes.
% INPUTS:
%  xml - jave object that will store xml data structure
%  parent - parent DOM Element
%  s - Matlab data structure to save
%  TagName - name to be used in xml tags describing 's'
%  Pref - preferenced
% OUTPUT:
%  parent - modified 'parent'

% perform some conversions
if (ischar(s) && min(size(s))>1) % if 2D array of characters
    s=cellstr(s);                  % than convert to cell array
end
% if (strcmp(TagName, 'CONTENT'))
%   while (iscell(s) && length(s)==1), s = s{1}; end % unwrap cell arrays of length 1
% end
TagName  = varName2str(TagName);

%% == node is a 2D cell array ==
% convert to some other format prior to further processing
nDim = nnz(size(s)>1);  % is it a scalar, vector, 2D array, 3D cube, etc?
if (iscell(s) && nDim==2 && strcmpi(Pref.CellTable, 'Matlab'))
    s = var2str(s, Pref.PreserveSpace);
end
if (nDim==2 && (iscell  (s) && strcmpi(Pref.CellTable,   'Vector')) || ...
        (isstruct(s) && strcmpi(Pref.StructTable, 'Vector')))
    s = s(:);
end
if (nDim>2), s = s(:); end % can not handle this case well
nItem = numel(s);
nDim  = nnz(size(s)>1);  % is it a scalar, vector, 2D array, 3D cube, etc?

%% == node is a cell ==
if (iscell(s)) % if this is a cell or cell array
    if ((nDim==2 && strcmpi(Pref.CellTable,'Html')) || (nDim< 2 && Pref.CellItem))
        % if 2D array of cells than can use HTML-like notation or if 1D array
        % than can use item notation
        if (strcmp(TagName, 'CONTENT')) % CONTENT nodes already have <TagName> ... </TagName>
            array2DOMnode(xml, parent, s, Pref.ItemName, Pref ); % recursive call
        else
            node = xml.createElement(TagName);   % <TagName> ... </TagName>
            array2DOMnode(xml, node, s, Pref.ItemName, Pref ); % recursive call
            parent.appendChild(node);
        end
    else % use  <TagName>...<\TagName> <TagName>...<\TagName> notation
        array2DOMnode(xml, parent, s, TagName, Pref ); % recursive call
    end
    %% == node is a struct ==
elseif (isstruct(s))  % if struct than deal with each field separatly
    if ((nDim==2 && strcmpi(Pref.StructTable,'Html')) || (nItem>1 && Pref.StructItem))
        % if 2D array of structs than can use HTML-like notation or
        % if 1D array of structs than can use 'items' notation
        node = xml.createElement(TagName);
        array2DOMnode(xml, node, s, Pref.ItemName, Pref ); % recursive call
        parent.appendChild(node);
    elseif (nItem>1) % use  <TagName>...<\TagName> <TagName>...<\TagName> notation
        array2DOMnode(xml, parent, s, TagName, Pref ); % recursive call
    else % otherwise save each struct separatelly
        fields = fieldnames(s);
        node = xml.createElement(TagName);
        for i=1:length(fields) % add field by field to the node
            field = fields{i};
            x = s.(field);
            switch field
                case {'COMMENT', 'CDATA_SECTION', 'PROCESSING_INSTRUCTION'}
                    if iscellstr(x)  % cell array of strings -> add them one by one
                        array2DOMnode(xml, node, x(:), field, Pref ); % recursive call will modify 'node'
                    elseif ischar(x) % single string -> add it
                        struct2DOMnode(xml, node, x, field, Pref ); % recursive call will modify 'node'
                    else % not a string - Ignore
                        warning('xml_io_tools:write:badSpecialNode', ...
                            ['Struct field named ',field,' encountered which was not a string. Ignoring.']);
                    end
                case 'ATTRIBUTE' % set attributes of the node
                    if (isempty(x)), continue; end
                    if (isstruct(x))
                        attName = fieldnames(x);       % get names of all the attributes
                        for k=1:length(attName)        % attach them to the node
                            att = xml.createAttribute(varName2str(attName(k)));
                            att.setValue(var2str(x.(attName{k}),Pref.PreserveSpace));
                            node.setAttributeNode(att);
                        end
                    else
                        warning('xml_io_tools:write:badAttribute', ...
                            'Struct field named ATTRIBUTE encountered which was not a struct. Ignoring.');
                    end
                otherwise                            % set children of the node
                    struct2DOMnode(xml, node, x, field, Pref ); % recursive call will modify 'node'
            end
        end  % end for i=1:nFields
        parent.appendChild(node);
    end
    %% == node is a leaf node ==
else  % if not a struct and not a cell than it is a leaf node
    switch TagName % different processing depending on desired type of the node
        case 'COMMENT'   % create comment node
            com = xml.createComment([s]);
            parent.appendChild(com);
        case 'CDATA_SECTION' % create CDATA Section
            cdt = xml.createCDATASection(s);
            parent.appendChild(cdt);
        case 'PROCESSING_INSTRUCTION' % set attributes of the node
            OK = false;
            if (ischar(s))
                n = strfind(s, ' ');
                if (~isempty(n))
                    proc = xml.createProcessingInstruction(s(1:(n(1)-1)),s((n(1)+1):end));
                    parent.insertBefore(proc, parent.getFirstChild());
                    OK = true;
                end
            end
            if (~OK)
                warning('xml_io_tools:write:badProcInst', ...
                    ['Struct field named PROCESSING_INSTRUCTION need to be',...
                    ' a string, for example: xml-stylesheet type="text/css" ', ...
                    'href="myStyleSheet.css". Ignoring.']);
            end
        case 'CONTENT' % this is text part of already existing node
            txt  = xml.createTextNode(var2str(s, Pref.PreserveSpace)); % convert to text
            parent.appendChild(txt);
        otherwise      % I guess it is a regular text leaf node
            txt  = xml.createTextNode(var2str(s, Pref.PreserveSpace));
            node = xml.createElement(TagName);
            node.appendChild(txt);
            parent.appendChild(node);
    end
end % of struct2DOMnode function

%% =======================================================================
%  === array2DOMnode Function ============================================
%  =======================================================================
function [] = array2DOMnode(xml, parent, s, TagName, Pref)
% Deal with 1D and 2D arrays of cell or struct. Will modify 'parent'.
nDim = nnz(size(s)>1);  % is it a scalar, vector, 2D array, 3D cube, etc?
switch nDim
    case 2 % 2D array
        for r=1:size(s,1)
            subnode = xml.createElement(Pref.TableName{1});
            for c=1:size(s,2)
                v = s(r,c);
                if iscell(v), v = v{1}; end
                struct2DOMnode(xml, subnode, v, Pref.TableName{2}, Pref ); % recursive call
            end
            parent.appendChild(subnode);
        end
    case 1 %1D array
        for iItem=1:numel(s)
            v = s(iItem);
            if iscell(v), v = v{1}; end
            struct2DOMnode(xml, parent, v, TagName, Pref ); % recursive call
        end
    case 0 % scalar -> this case should never be called
        if ~isempty(s)
            if iscell(s), s = s{1}; end
            struct2DOMnode(xml, parent, s, TagName, Pref );
        end
end

%% =======================================================================
%  === var2str Function ==================================================
%  =======================================================================
function str = var2str(object, PreserveSpace)
% convert matlab variables to a string
switch (1)
    case isempty(object)
        str = '';
    case (isnumeric(object) || islogical(object))
        if ndims(object)>2, object=object(:); end  % can't handle arrays with dimention > 2
        str=mat2str(object);           % convert matrix to a string
        % mark logical scalars with [] (logical arrays already have them) so the xml_read
        % recognizes them as MATLAB objects instead of strings. Same with sparse
        % matrices
        if ((islogical(object) && isscalar(object)) || issparse(object)),
            str = ['[' str ']'];
        end
        if (isinteger(object)),
            str = ['[', class(object), '(', str ')]'];
        end
    case iscell(object)
        if ndims(object)>2, object=object(:); end  % can't handle cell arrays with dimention > 2
        [nr nc] = size(object);
        obj2 = object;
        for i=1:length(object(:))
            str = var2str(object{i}, PreserveSpace);
            if (ischar(object{i})), object{i} = ['''' object{i} '''']; else object{i}=str; end
            obj2{i} = [object{i} ','];
        end
        for r = 1:nr, obj2{r,nc} = [object{r,nc} ';']; end
        obj2 = obj2.';
        str = ['{' obj2{:} '}'];
    case isstruct(object)
        str='';
        warning('xml_io_tools:write:var2str', ...
            'Struct was encountered where string was expected. Ignoring.');
    case isa(object, 'function_handle')
        str = ['[@' char(object) ']'];
    case ischar(object)
        str = object;
    otherwise
        str = char(object);
end

%% string clean-up
str=str(:); str=str.';            % make sure this is a row vector of char's
if (~isempty(str))
    str(str<32|str==127)=' ';       % convert no-printable characters to spaces
    if (~PreserveSpace)
        str = strtrim(str);             % remove spaces from begining and the end
        str = regexprep(str,'\s+',' '); % remove multiple spaces
    end
end

%% =======================================================================
%  === var2Namestr Function ==============================================
%  =======================================================================
function str = varName2str(str)
% convert matlab variable names to a sting
str = char(str);
p   = strfind(str,'0x');
if (~isempty(p))
    for i=1:length(p)
        before = str( p(i)+(0:3) );          % string to replace
        after  = char(hex2dec(before(3:4))); % string to replace with
        str = regexprep(str,before,after, 'once', 'ignorecase');
        p=p-3; % since 4 characters were replaced with one - compensate
    end
end
str = regexprep(str,'_COLON_',':', 'once', 'ignorecase');
str = regexprep(str,'_DASH_' ,'-', 'once', 'ignorecase');


% -------------------load_sto_file ----------------------------------%
function out = load_sto_file(filename)

% function data = load_sto_file(filename,delimiters)
%
% This function loads either STO or MOT files and stores each column of
% data in structure array with the field name which is the columning
% heading in the file. It discards any other information from the header of 
% the file.
%
% Input: filename - the STO or MOT filename
%
% Author: Glen Lichtwark 
% Last Modified: 17/11/2008
%   Updates:
%   L71 -  ~isempty(str2num(f_name(1)))&& ~contains(f_name,'iliacus') to
%   avoid the iliacus to be replaced by Niliacus (BG, 2/10/2020)

if nargin < 1
    [fname, pname] = uigetfile('*.*', 'File to load - ');
    file = [pname fname];
else file = filename;    
end

[file_data,s_data]= readtext(file, '\t', '', '', 'empty2NaN');

% search the numerical data (in s_data.numberMask) to find when the block 
% of data starts

a = find(abs(diff(sum(s_data.numberMask,2)))>0);
[m,n] = size(file_data);

% create an array with all of the data
num_dat = [file_data{a(end)+1:end,1:sum(s_data.numberMask(a(end)+1,:),2)}];

% reshape to put back into columns and rows
data = reshape(num_dat,m-a(end),sum(s_data.numberMask(a(end)+1,:),2));

% now find the column headings if there are any
if sum(s_data.stringMask(a(end),:)) == sum(s_data.numberMask(a(end)+1,:))
    data_label = file_data(a(end),:);
    b = a(end)-1;
else
    cols = sum(s_data.numberMask(a(end)+1,:));
    data_label = file_data(a(end),1:cols);
    b = a(end);
end

% go through the data labels and find any that are duplicates (this occurs
% in the ground reaction force data where each forceplate has the same
% column headings) and add a number to distinguish the duplicates.

for i = 1:length(data_label)
    tf = strcmp(data_label(i),data_label);
    c = find(tf>0);
    if length(c) > 1
        for j = 1:length(c)
            data_label(c(j)) = cellstr([data_label{c(j)} num2str(j)]);
        end
    end
end

% now create the output structure with the field names from the data labels
% and the corresponding data from the columns of the data array
for i = 1:length(data_label)
    f_name = data_label{i};
    % find any spaces and replace with underscore
    e = findstr(f_name, ' ');
    if ~isempty(e)
        f_name(e) = '_';
    end
    e = findstr(f_name, '.');
    if ~isempty(e)
        f_name(e) = '_';
    end
    if ~isempty(str2num(f_name(1)))&& ~contains(f_name,'iliacus')
        f_name = ['N' f_name];
    end
    out.(f_name) = data(:,i);
end


% ----------------------- readtext ---------------------------%
function      [data, result]= readtext(fname, delimiter, comment, quotes, options)

%  Usage:     [data, result]= readtext(fname, delimiter, comment, quotes, options)
% 
% Whatever text file you give it, readtext returns an array of the contents (or send me a 
%   bug report). Matlab can't read variable length lines or variable type values with the standard 
%   library. readtext can read any text file. Any string (or even regexp) can be delimiting, 
%   default is a comma. Everything after (and including) a comment character, until the line end, 
%   is ignored. Quote characters may also be given, everything between them is treated as one item. 
%   There are options to control what will be converted to numbers and how empty items are saved. 
% 
% If you find any errors, please let me know: peder at axensten dot se
% 
% fname:      the file to be read.
% 
% delimiter:  (default: ',') any string. May be a regexp, but this is a bit slow on large files. 
% 
% comment:    (default: '') zero or one character. Anything after (and including) this character, 
%   until the end of the line, will be ignored. 
% 
% quotes:     (default: '') zero, one (opening quote equals closing), or two characters (opening 
%   and closing quote) to be treated as paired braces. Everything between the quotes will be 
%   treated as one item. The quotes will remain. Quotes may be nested.
% 
% options:    (default: '') may contain (concatenate combined options): 
% - 'textual': no numeric conversion ('data' is a cell array of strings only), 
% - 'numeric': everything is converted to a number or NaN ('data' is a numeric array, empty items 
%   are converted to NaNs unless 'empty2zero' is given), 
% - 'empty2zero': an empty field is saved as zero, and 
% - 'empty2NaN': an empty field is saved as NaN. 
% - 'usewaitbar': call waitbar to report progress. If you find the wait bar annoying, get 'waitbar 
%   alternative' at http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=11398
% 
% data:       A cell array containing the read text, divided into cells by delimiter and line 
%   endings. 'data' will be empty if the file is not found, could not be opened, or is empty. 
%   With the option 'numeric', 'data' will be a numeric array, with 'textual', 'data' will be a 
%   cell array of strings only, and otherwise it will be a mixed cell array. For Matlab < version 7, 
%   returned strings may contain leading white-space.
% 
% result:     a structure:
% .min: minimum number of columns found in a line.
% .max: number of columns in 'data', before removing empty columns.
% .rows: number of rows in 'data', before removing empty rows. 
% .numberMask: true, if numeric conversion ('NaN' converted to NaN counts).
% .number: number of numeric conversions ('NaN' converted to NaN counts).
% .emptyMask: true, if empty item in file.
% .empty: number of empty items in file.
% .stringMask: true, if non-number and non-empty.
% .string: number of non-number, non-empty items.
% .quote: number of quotes. 
% 
% INSPIRATION: loadcell.m (id 1965). The goal of readtext is to be at least as flexible (you be 
%   the judge) and quicker. On my test file (see below), readtext is about 3--4 times 
%   as quick, maybe even more on large files. In readtext you may use a regexp as 
%   delimiter and it can ignore comments in the text file. 
% 
% SPEED:      Reading a 1MB file (150000 items!) with 'numeric' takes about 100 seconds on a 
%   fairly slow system. Time scales approximately linearly with input file size. 
% - Conversion from string to numeric is slow (I can't do anything about this), but using the 
%   option 'textual' is a lot quicker (above case takes 12 seconds).
% - Using a regexp delimiter is slower (during initializing), it adds 250 seconds! 
% 
% EXAMPLE:    [a,b]= readtext('txtfile', '[,\t]', '#', '"', 'numeric-empty2zero')
% This will load the file 'txtfile' into variable a, treating any of tab or comma as
%   delimiters. Everything from and including # to the next newline will be ignored. 
%   Everything between two double quotes will be treated as a string. Everything will 
%   be converted to numbers and a numeric array returned. Non-numeric items will become 
%   NaNs and empty items are converted to zero. 
% 
% Copyright (C) Peder Axensten (peder at axensten dot se), 2006.

% HISTORY:
% Version 1.0, 2006-05-03.
% Version 1.1, 2006-05-07:
% - Made 'options' case independent. 
% - Removed result.errmess -- now use error(errmess) instead. 
% - Removed result.nan -- it was equivalent to result.string, so check this instead.
% - Added some rows', 'result' fields: 'numberMask', 'emptyMask', and 'stringMask' 
%   (see 'result:', above).
% - A few small bug fixes.
% Version 1.2, 2006-06-06:
% - Now works in Matlab 6.5.1 (R13SP1) (maybe 6.5 too), versions <6.5 will NOT work.
% Version 1.3, 2006-06-20:
% - Better error report when file open fails. 
% - Somewhat quicker. 
% - Recommends 'waitbar alternative'. Ok with Matlab orig. waitbar too, of course. 
% Version 1.4, 2006-07-14:
% - Close waitbar instead of deleting it, and some other minor waitbar compatibility fixes. 
% Version 1.5, 2006-08-13:
% - No more (La)TeX formatting of file names. 
% - Prefixed waitbar messages with '(readtext)'. 
% Version 1.6, 2006-10-02:
% - Better removal of comments. Could leave an empty first row before. 
% - Added a 'usewaitbar' option. 
% - Now removes empty last columns and rows. 
% 
% TO DO:
% - Write a better text macro. 
% - Add result.quoteMask.
% - Add 'removeemptycolumns' and 'removeemptyrows' options. 

% KEYWORDS:     import, read, load, text, delimited, cell, numeric, array, flexible
% 
% REQUIREMENTS: Works in Matlab 6.5.1 (R13SP1) (probably 6.5 too), versions <6.5 will NOT work.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Read (or set to default) the input arguments.
	if((nargin < 1) || ~ischar(fname) || isempty(fname))		% Is there a file name?
		error('First argument must be a file name!'); 
	end
	if(nargin < 2), delimiter=	',';				end			% Default delimiter value.
	if(nargin < 3), comment=	'';					end			% Default comment value.
	if(nargin < 4), quotes=		'';					end			% Default quotes value.
	if(nargin < 5), options=	[];					end			% Default options value.
	
	options=		lower(options);
	op_waitbar=		~isempty(strfind(options, 'usewaitbar'));	% Do waitbar calls. 
	op_numeric=		~isempty(strfind(options, 'numeric'));		% Read as numerical. 
	op_textual=		~isempty(strfind(options, 'textual')) && ~op_numeric;	% Read as textual. 
	op_empty=		[];											% Ignore empties, ...
	if(~isempty(strfind(options, 'empty2zero')))
		op_empty=		0;										% ... or replace by zero ...
	elseif(op_numeric || ~isempty(strfind(options, 'empty2nan')))
		op_empty=		NaN;									% ... or replace by NaN.
	end
	if(op_textual), op_empty= num2str(op_empty);	end			% Textual 'empty'.
	if(~ischar(comment) || (length(comment) > 1))
		error('Argument ''comment'' must be a string of maximum one character.');
	end
	if(~ischar(quotes) || (length(quotes) > 2))
		error('Argument ''quotes'' must be a string of maximum two characters.');
	end
	
	% Set the default return values.
	result.min=		Inf;
	result.max=		0;
	result.quote=	0;
	
	% Read the file.
	[fid, errmess]=	fopen(fname, 'r');							% Open the file.
	if(fid < 0), error(['Trying to open ' fname ': ' errmess]); end
	text=			fread(fid, 'uchar=>char')';					% Read the file.
	fclose(fid);												% Close the file.
	
	if(op_waitbar)
		th= waitbar(0, '(readtext) Initialising...');			% Show waitbar.
		thch=			findall(th, '-property', 'Interpreter');
		set(thch, 'Interpreter', 'none');						% No (La)TeX) formatting. 
	end
	
	% Clean up the text.
	eol=			char(10);
	text=			strrep(text, [char(13) char(10)], eol);		% Replace Windows-style eol.
	text=			strrep(text, char(13), eol);				% Replace MacClassic-style eol.
	if(~isempty(comment))										% Remove comments.
		text=	regexprep(text, ['^' comment '[^\n]*\n'], '');	% Remove commented lines. 
		text=	regexprep(text, [comment '[^\n]*'], '');		% Remove commented line endings. 
	end
	if(text(end) ~= eol), text= [text eol];				end		% End string with eol, if none.
	
	% Find column and row dividers.
	delimiter=		strrep(delimiter, '\t', char( 9));			% Convert to one char, quicker?
	delimiter=		strrep(delimiter, '\n', char(10));
	delimiter=		strrep(delimiter, '\r', char(10));
	delimiter=		strrep(delimiter, '\f', char(12));
	if(1 == length(delimiter))									% Find column dividers quickly.
		delimS=		find((text == delimiter) | (text == eol));
		delimE=		delimS;
	elseif(isempty(regexp(delimiter, '[\+\*\?\|\[^$<>]', 'once')))	% Find them rather quickly.
		delimS=		strfind(text, delimiter);
		eols=		find(text == eol);
		delimE=		union(eols, delimS + length(delimiter) - 1);
		delimS=		union(eols, delimS);
	else														% Find them with regexp.
		[delimS, delimE]=	regexp(text, [delimiter '|' eol]);
	end
	divRow=			[0, find(text == eol), length(text)];		% Find row dividers+last.
	
	% Keep quoted text together.
	if(~isempty(quotes))										% Should we look for quotes?
		if((length(quotes) == 1) || (quotes(1) == quotes(2)))	% Opening char == ending.
			exclE=			find(text == quotes(1));
			exclS=			exclE(1:2:end);
			exclE=			exclE(2:2:end);
		else													% Opening char ~= closing.
			exclS=			find(text == quotes(1));
			exclE=			find(text == quotes(2));
		end
		if((length(exclS) ~= length(exclE)) || (sum(exclS > exclE) > 0))
			if(op_waitbar), close(th); 	end						% Close waitbar or it'll linger.
			error('Opening and closing quotes don''t match in file %s.', fname); 
		end
		if(~isempty(exclS))										% We do have quoted text.
			if(op_waitbar), waitbar(0, th, '(readtext) Doing quotes...'); end	% Inform user.
			r=		1;
			rEnd=	length(exclS);
			n=		1;
			nEnd=	length(delimS);
			result.quote=	rEnd;
			while((n < nEnd) && (r < rEnd)) % "Remove" delimiters and newlines within quyotes.
				while((r <= rEnd) && (delimS(n) > exclE(r))), r= r+1;	end
				while((n <= nEnd) && (delimS(n) < exclS(r))), n= n+1;	end
				while((n <= nEnd) && (delimS(n) >= exclS(r)) && (delimS(n) <= exclE(r)))
					delimS(n)=	0;
					n=			n+1;
				end
				if((bitand(n, 255) == 0) && op_waitbar), waitbar(n/nEnd); end	% Update waitbar.
			end
			if(op_waitbar), waitbar(1);	end;
			delimE=	delimE(delimS > 0);
			delimS=	delimS(delimS > 0);
		end
	end
	delimS=		delimS-1;										% Last char before delimiter.
	delimE=		[1 delimE(1:end-1)+1];							% First char after delimiter.
	
	% Do the stuff: convert text to cell (and maybe numeric) array.
	if(op_waitbar), waitbar(0, th, sprintf('(readtext) Reading ''%s''...', fname));	end
	r=				1;
	c=				1;											% Presize data to optimise speed.
	data=			cell(length(divRow), ceil(length(delimS)/(length(divRow)-1)));
	nums=			zeros(size(data));							% Presize nums to optimise speed.
	nEnd=			length(delimS);								% Prepare for a waitbar.
	for n=1:nEnd
		temp=			text(delimE(n):delimS(n));
		data{r, c}= 	temp;									% Textual item.
		if(~op_textual), nums(r, c)= str2double(temp);	end		% Quicker(!) AND better waitbar.
		if(text(delimS(n)+1) == eol)							% Next row.
			result.min=		min(result.min, c);					% Find shortest row.
			result.max=		max(result.max, c);					% Find longest row.
			r=				r+1;
			c=				0;
		end
		c=				c+1;
		if((bitand(n, 255) == 0) && op_waitbar), waitbar(n/nEnd);	end	% Update waitbar.
	end
	if(op_waitbar), waitbar(1);	end
	
	% Clean up the conversion and do the result statistics.
	if(op_waitbar), waitbar(0, th, '(readtext) Cleaning up...');	end	% Inform user.
	data=				data(1:(r-1), 1:result.max);			% In case we started off to big.
	if(~op_textual), nums= nums(1:(r-1), 1:result.max);	end		% In case we started off to big.
	if(exist('strtrim', 'builtin')), data= strtrim(data);		% Not in Matlab 6.5...
	else							 data= deblank(data);		
	end
	while(all(cellfun('isempty', data(end, :))))				% Remove empty last lines. 
		data=	data(1:end-1, :); 
		nums=	nums(1:end-1, :); 
		r=		r-1;
	end 
	while(all(cellfun('isempty', data(:, end))))				% Remove empty last columns. 
		data=	data(:, 1:end-1); 
		nums=	nums(:, 1:end-1); 
		c=		c-1;
	end 
	result.rows=		r-1;
	empties=			cellfun('isempty', data);				% Find empty items.
	result.emptyMask=	empties;
	if(op_textual)
		result.numberMask=	repmat(false, size(data));			% No numbers, all strings.
		result.stringMask=	~empties;							% No numbers, all strings.
		data(empties)=		{op_empty};							% Set correct empty value.
	else
		result.numberMask=	~(isnan(nums) & ~strcmp(data, 'NaN'));	% What converted well.
		if(op_numeric)
			nums(empties)=		op_empty;						% Set correct empty value.
			data=				nums;							% Return the numeric array.
			result.stringMask=	~(empties | result.numberMask);	% Didn't convert well: so strs.
		else
			data(result.numberMask)= num2cell(nums(result.numberMask));	% Copy back numerics.
			data(empties)=		{op_empty};						% Set correct empty value.
			result.stringMask=	cellfun('isclass', data, 'char');	% Well, the strings.
		end
	end
	result.empty=		sum(result.emptyMask(:));				% Count empties.
	result.numberMask=	result.numberMask & ~result.emptyMask;	% Empties don't count.
	result.number=		sum(result.numberMask(:));				% Count numbers.
	result.stringMask=	result.stringMask & ~result.emptyMask;	% Empties don't count.
	result.string=		sum(result.stringMask(:));				% Count strings.
	
 	if(op_waitbar), close(th);	end								% Removing the waitbar. 





