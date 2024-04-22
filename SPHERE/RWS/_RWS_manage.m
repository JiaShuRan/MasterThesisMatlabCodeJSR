
function RWS_manage

files = dir('sf*');

for id = 1:length(files)
filesL=files(id);
filename=filesL.name;
filename5=strcat(extractBefore(filename,6),'.dat');
movefile(filename,filename5)

end




