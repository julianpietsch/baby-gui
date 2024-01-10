function deleteFileAnnotations(obj, fileAnnotations)
%Deletes the input file annotation (including the original file and the
%link)
hm=java.util.HashMap;
for n=1:length(fileAnnotations)
    l=java.util.LinkedList;
    l.add(java.lang.Long(fileAnnotations(n).getId.getValue));
    del=omero.cmd.Delete2;
    hm.put('Annotation',l);
    del.targetObjects=hm;
    obj.Session.submit(del);
end

end