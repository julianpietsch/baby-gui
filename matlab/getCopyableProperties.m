function pnames = getCopyableProperties(obj,classname)
%GETCOPYABLEPROPERTIES Get only copyable properties for an object
%   PNAMES = getCopyableProperties(OBJ,CLASSNAME)
%   If provided with an object of class CLASSNAME (or one of its
%   subclasses), this function returns all properties, both public and
%   private, that can be copied by methods in CLASSNAME classes, but omits
%   properties that should not be set (i.e., Transient, Dependent, Constant
%   and immutable properties). For structures, it returns all of the
%   fieldnames.
%   
%   This basic helper function gets used in the "copyprops" method of
%   babyExperiment and babyTimelapse classes and their children,  which in
%   turn is called by the "loadobj" method for each of those classes.
%
%   See also BABYEXPERIMENT and BABYTIMELAPSE

if ~ischar(classname) || ~isvector(classname)
    error('"classname" must be a char vector');
end

if isstruct(obj)
    pnames = fieldnames(obj);
else
    % Get all properties accessible by this class
    m = meta.class.fromName(classname);
    if isa(obj,classname)
        props = m.PropertyList;
    elseif length(m.SuperclassList)==1 && isa(obj,m.SuperclassList.Name) && ...
            ~strcmp(m.SuperclassList.Name,'handle')
        props = m.SuperclassList.PropertyList;
    else
        % The following should never happen, but will get raised as a warning
        % if it is encountered in "loadobj" methods:
        error(['This does not look like a "' classname '" object.']);
    end
    setAccess = {props.SetAccess};
    validSetAccess = cellfun(@ischar,setAccess);
    setAccess(~validSetAccess) = {'complicated'};
    pnames = {props.Name};
    pnames = unique(pnames(~[props.Transient] & ...
        ~[props.Dependent] & ~[props.Constant] & ...
        ~ismember(setAccess,{'immutable','none'})));
end

end
