%Helper file for Belham_Valley_datareader

function [yada, xyz] =read_in_zgrid2(filename)
%filename
fid=fopen(filename,'r');
%fid
yada=fscanf(fid,'Nx=%g: X={%g,%g}\nNy=%g: Y={%g,%g',[3 2])
fgets(fid); fgets(fid); %'}\nElevation=\n' OR '}\nPileheight=\n'

xyz=zeros(yada(1,1),yada(1,2),3);
size(xyz);

xyz(:,:,1)=((2*(0:(yada(1,1)-1))+0.5)/(2*yada(1,1))*(yada(3,1)-yada(2,1))+yada(2,1))'*ones(1,yada(1,2)); %'; get one column and repmat it into yada(1,2) columns 


xyz(:,:,2)=ones(yada(1,1),1)*(2*(0:(yada(1,2)-1))+0.5)/(2*yada(1,2))*(yada(3,2)-yada(2,2))+yada(2,2); % get one row a repmat it into a repmat it into yada(1, 1) rows



xyz(:,:,3)=fscanf(fid,'%g',[yada(1,1) yada(1,2)]);
fclose(fid);

end

