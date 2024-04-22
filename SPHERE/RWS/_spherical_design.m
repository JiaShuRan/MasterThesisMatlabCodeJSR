
function xyzw=spherical_design(deg)

%--------------------------------------------------------------------------
% Object:
%--------------------------------------------------------------------------
% This routine loads spherical designs computed by Rob Womersley, useful
% for cubature on the unit 2-sphere, with algebraic degree of precision at 
% most "deg", with "deg <= 180".
%--------------------------------------------------------------------------
% Input:
%--------------------------------------------------------------------------
% deg: algebraic degree of precision of the cubature rule.
%--------------------------------------------------------------------------
% Output:
%--------------------------------------------------------------------------
% xyzw: nodes and weights are stored in a "M x 4" matrix, where "M" is the
%     cardinality of the rule.
%     The k-th node, in cartesian coordinates correspond to "xyzw(k,1:3)",
%     while the k-th weight is "xyzw(k,4)".
%     The weights are choosen so that their sum is "4*pi"
%--------------------------------------------------------------------------
% Reference:
%--------------------------------------------------------------------------
% All these rules are taken from 
% "Efficient Spherical Designs with Good Geometric Properties"
% and are studied and stored by Rob Womersley. 
%--------------------------------------------------------------------------

if deg >= 180
    fprintf(2,'\n \t The ADE is too high. Setting deg=180.');
    deg=180;
end

degstr=num2str(deg);
if deg < 100
    degstr=strcat('0',degstr);
end

filename=strcat('sf',degstr,'.dat');
files=dir(filename)

fName=files.name;
fFolder=files.folder;
fullname=fullfile(fFolder,fName)

xyzwVALS=load(fullname);

xyzw=xyzwVALS;