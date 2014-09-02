function [subDirs]=List_Folders(url_catalog)
% List_Folders needs an http access to download the xml file of the
% thredds catalog. This function lists all the subFolders of a Thredds catalog page.
%
% example : 
% url_catalog='http://opendap-qcif.arcs.org.au/thredds/catalog/IMOS/ACORN/gridded_1h-avg-current-map_non-QC/CBG/2007/catalog.xml';
%
% Inputs:
%   url_catalog                - https address of the THREDDS catalog
%                                .XML,Not HTML
% 
% Outputs :
%   subdir                   - Cell array of the folders
%
% Author: Laurent Besnard <laurent.besnard@utas.edu.au>
%
%
% Copyright (c) 2011, eMarine Information Infrastructure (eMII) and Integrated
% Marine Observing System (IMOS).
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%     * Redistributions of source code must retain the above copyright notice,
%       this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the eMII/IMOS nor the names of its contributors
%       may be used to endorse or promote products derived from this software
%       without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%

warning off all
%% Put in a structure called V the content of url_catalog
url_catalog={url_catalog};
fileList=[]';

filenameXML='THREDDS.xml';
urlwrite(url_catalog{1}, filenameXML);
V = xml_parseany(fileread(filenameXML));
delete(filenameXML);


%% List the subfolders in the current page
try
    NumberSubfolders=length(V.dataset{1}.catalogRef);
    subDirs=cell(NumberSubfolders,1);
    for ii=1:NumberSubfolders
        subDirs{ii}=V.dataset{1}.catalogRef{ii}.ATTRIBUTE.xlink_href;
        
        
    end
catch
    subDirs=[];
    NumberSubfolders=0;
end