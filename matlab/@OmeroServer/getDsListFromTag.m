function [datasets,dsNames]=getDsListFromTag(this, tagText)
%Returns arrays of dataset names and objects tagged with the input text
tags = getObjectAnnotations(this.session,'tag','dataset',[]);
tagnames = arrayfun(@(t) char(t.getTextValue.getValue),tags,'uni',0);
tagids = arrayfun(@(t) t.getId.getValue,tags);
tagids = unique(tagids(contains(tagnames,tagText,'IgnoreCase',true)));

if isempty(tagids)
    datasets = [];
    dsNames = {};
    return
end

% Following is based on omero-py/src/omero/gateway/__init__.py line 4174
% in the getObjectsByAnnotations function
sql = [...
    'select ob from Dataset ob ',...
    'left outer join fetch ob.annotationLinks obal ',...
    'left outer join fetch obal.child ann ',...
    'where ann.id in (:oids)'];

parameters = omero.sys.ParametersI();
oids = omero.rtypes.rlist(toJavaList(tagids,@omero.rtypes.rlong));
parameters.map.put('oids', oids);
context = java.util.HashMap;

proxy = this.session.getQueryService();
dsList = proxy.findAllByQuery(sql,parameters,context);
datasets = toMatlabList(dsList);
dsNames = arrayfun(@(d) char(d.getName.getValue),datasets,'uni',0);
end
