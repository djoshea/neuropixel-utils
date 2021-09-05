function out = strings2java(str)
out = javaArray('java.lang.String', numel(str));
for i = 1:numel(str)
  out(i) = java.lang.String(str(i));
end
end