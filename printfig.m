function printfig(h, filename, varargin)
%printfig: Funktion zum Export von Figures als Vektor- oder
%Bitmap-Grafiken. Die Paperpostion wird angepasst, sodass Vektorgrafiken
%moeglichst korrekt abgespeichert werden.
%
%Syntax:
%
%       printfig(h, 'filename')
%       printfig(h, 'filename', 'filetype1', 'filetype2')
%       printfig(h, 'filename', 'filetype1, option1a, option1b, 'filetype2')
%       printfig(h, 'filename', 'expandaxes', 'filetype')
%       printfig(h, 'filename', 'expandaxes', fHor, fVer, 'filetype1', 'filetype2')
%       printfig(h) %zum Kopieren (Clipboard)
%
%Bsp.:  printfig(gcf, 'Beispiel', 'emf', 'png', 1200, 'transparent', 'expandaxes', 1.5, 1.5)
%
%Eingabeparameter:
%
%   h:          Figure-Handle
%   filename:   Dateiname als String (Ohne Dateiendung)
%
%
%
%Optionale Inputs:
%
%   Datentyp (mehrere gleichzeitig moeglich) zum speichern:
%           'emf','eps','bmp','jpg','jpeg','tiff','tif','pdf','png','copy', 'svg' und/oder 'fig'
%           Wird kein Datentyp vorgegeben, wird das Figure als emf gespeichert.
%           Bei Eingabe von 'copy' wird eine Vektorgrafik in die
%           Zwischenablage kopiert
%           Wird kein Dateiname vorgegeben, wird das Figure nur in die Zwischenablage
%           kopiert.
%
%   DPI (hinter den Datentyp) - default: Bildschirmaufl�sung (i.d.R. 600 dpi)
%       z. B. printfig(h,'BspName', 'png', 1200, 'jpg', 600, 'emf')
%
%   'tpb'           Einzelne Bilder werden mit transparentem Hintergrund
%                   gespeichert, falls moeglich. Diese Option kann hinter
%                   jede Datei gesetzt werden.
%
%                   Bsp.: printfig(gcf, 'tiff', 1200, 'png', 'tpb')
%                   --> speichert eine TIFF-Datei mit 1200 DPI und nicht
%                   transparentem Hintergrund und eine PNG-Datei mit
%                   transparentem Hintergrund
%
%   'transparent'   Alle Bilder werden mit transparentem Hintergrund
%                   gespeichert, falls moeglich. 'tpb'-Inputs werden in
%                   diesem Fall ignoriert.
%
%   'expandaxes'    Es wird versucht, expandaxes() auf das Figure anzuwenden
%                   Hinter die Option 'expandaxes' koennen die
%                   Abstandsfaktoren fHor und fVer vorgegeben werden (doc
%                   expandaxes fuer Details)
%                   Bsp.: printfig(gcf, 'expandaxes', 'png')
%                         printfig(gcf, 'expandaxes', 1.5, 1.5, 'png')
%
%
%
%Ausgabeparameter:
%
%   keine
%
%
%verwendete Funktionen:
%
%   expandaxes.m (wenn 'expandaxes'-Option verwendet wird
%
%
%Autor: Marc Jakobi 08.07.2016
%
%
%aktuell bekannte Bugs:
%
%                     - Bitmap-Dateien und das Clipboard-Kopieren sind mit
%                       transparenten Optionen noch nicht moeglich.
%                       Transparente Hintergruende sind derzeit nur mit
%                       Vektorgrafiken moeglich.

%% Changelog:
%   13.07.2016:       - Option hinzugefuegt, Grafiken mit transparentem
%                       Hintergrund zu speichern.
%                     - Option hinzugefuegt, fHor und fVer
%                       (expandaxes-Parameter) vorzugeben.
%                     - Bei mehrfach gleichen Inputs wird nur die erste
%                       Konfiguration verwendet.
%   22.07.2016        - Moeglichkeit hinzugefuegt, Vektorgrafiken mit
%                       transparentem Hintergrund zu speichern (bitmap
%                       funktioniert noch nicht mit transparentem
%                       Hintergrund)
%   03.08.2016        - Option, als svg-Vektorgrafik zu speichern
%                       hinzugefuegt
%                     - Unterstuetzung fuer UNIX-Systeme hinzugefuegt:
%                       --> Clipboard-Copy als svg (oder pdf, wenn svg
%                           fehlschlaegt)
%                       --> svg fuer UNIX-Systeme als Default, statt emf

fHor = 1; fVer = 1;
%% Parse inputs
if ~isa(h,'matlab.ui.Figure')
    error('Der erste Input muss ein Figure-Handle sein')
end
if nargin == 1 %no file name
    varargin = {'copy'};
    filename = 'none';
elseif nargin == 2 
    if ispc %default: emf (Windows)
        varargin = {'emf'};
    else %default: svg (Unix)
        varargin = {'svg'};
    end
end

%check for file extension in file name and delete
if numel(filename) > 4 && strcmp(filename(end-3),'.')
    filename = filename(1:end-4);
end
%list of possible file types + expandaxes option
expectedFileTypes = {'expandaxes','emf','eps','bmp','jpg','jpeg','tiff','tif','pdf','png','copy','fig','svg'};
FTs = false(size(expectedFileTypes))'; %chosen file types
res = zeros(size(FTs)); %resolution vector (respective to chosen file types)
tpb = FTs; %transparent backgrounds
for i = 1:length(FTs)
    idx = StrfindInCell(varargin,expectedFileTypes{i},2);
    if ~isempty(idx)
        if length(idx) > 1 %more than one input of the same file type?
            warning(['Mehrere Inputs f�r ''',varargin{idx(1)},''' gefunden. ',...
                'Es wird nur die erste Konfiguration beachtet.'])
            varargin{idx(2:end)} = ' '; %delete later inputs
            idx = idx(1);
        end
        FTs(i) = ~isempty(idx);
        %non-default resolutions?
        if isnumeric(varargin{min(idx+1,end)})
            res(i) = varargin{min(idx+1,end)};
        end
        %transparent backgrounds?
        if strcmp(varargin{min(idx+1,end)},'tpb') || ...
                strcmp(varargin{min(idx+2,end)},'tpb')
            tpb(i) = true;
        end
        %check expandaxes option for fHor or fVer inputs
        if strcmp(varargin{idx},'expandaxes')
            if isnumeric(varargin{min(idx+1,end)})
                fHor = varargin{min(idx+1,end)};
                if isnumeric(varargin{min(idx+2,end)})
                    fVer = varargin{min(idx+1,end)};
                else
                    fVer = 1;
                end
            else
                fHor = 1;
            end
        end
    end
end
%if no type chosen
if sum(FTs) == 0
    FTs(1) = true;
elseif FTs(4) && FTs(5)
    warning('''jpg'' und ''jpeg'' wurden vorgegeben. Es kann sein, dass die ''jpeg''-Datei die ''jpg''-Datei �berschreibt.')
end
if FTs(6) && FTs(7)
    warning('''tif'' und ''tiff'' wurden vorgegeben. Es kann sein, dass die ''tiff''-Datei die ''tif''-Datei �berschreibt.')
end


%% Undock docked figures temporarily
WindowStyle = h.WindowStyle;
h.WindowStyle = 'normal';

%% Axes objects
AX = findobj(h,'type','axes');

%% transparent background?
%original colour values
bckcO = h.Color; %background color
colO = getaxbcol(AX); %sub-function to return axes background colours in struct
renderer = h.Renderer;
%if yes, set colour to none before saving colour to override other settings
if ~isempty(StrfindInCell(varargin,'transparent',2))
   h.Color = 'none';
   setaxbcol(AX,'none') %sub-function to set color of all elements of AX
   set(h,'InvertHardcopy','off')
   set(h,'renderer','painter')
end
bckc = h.Color; %background color
col = getaxbcol(AX); %sub-function to return axes background colours in struct

%% expandaxes?
if FTs(1) %expandaxes
    expandaxes(h, fHor, fVer)
end

%% Correct paper (printing) positions for vector graphics export
unis = get(h,'units');
ppos = get(h,'paperposition');
set(h,'units',get(h,'paperunits'));
pos = get(h,'position');
ppos(3:4) = pos(3:4);
set(h,'paperposition',ppos);
set(h,'units',unis);


%% save files
if FTs(2) %emf
    if tpb(2) %transparent background?
        h.Color = 'none';
        setaxbcol(AX,'none')
        set(h,'InvertHardcopy','off')
        set(h,'renderer','painter');
    end
    print(h,[filename,'.emf'],'-dmeta');
end
if FTs(3) %eps
    if tpb(3) %transparent background?
        h.Color = 'none';
        setaxbcol(AX,'none')
        set(h,'InvertHardcopy','off')
        set(h,'renderer','painter');
    else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
    end
    print(h,[filename,'.eps'],'-depsc')
end
if FTs(4) %bmp
%     if tpb(4) %transparent background?
%         h.Color = 'none';
%         setaxbcol(AX,'none')
%         set(h,'InvertHardcopy','off')
%         set(h,'renderer','painter');
%     else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
%     end
    print(h,[filename,'.bmp'],'-dbmp',['-r',num2str(res(4))]);
end
if FTs(5) %jpg
%     if tpb(5) %transparent background?
%         h.Color = 'none';
%         setaxbcol(AX,'none')
%         set(h,'InvertHardcopy','off')
%         set(h,'renderer','painter');
%     else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
%     end
    print(h,[filename,'.jpg'],'-djpeg',['-r',num2str(res(5))]);
end
if FTs(6) %jpeg
%     if tpb(6) %transparent background?
%         h.Color = 'none';
%         setaxbcol(AX,'none')
%         set(h,'InvertHardcopy','off')
%         set(h,'renderer','painter');
%     else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
%     end
    print(h,[filename,'.jpeg'],'-djpeg',['-r',num2str(res(6))]);
end
if FTs(7) %tif
%     if tpb(7) %transparent background?
%         h.Color = 'none';
%         setaxbcol(AX,'none')
%         set(h,'InvertHardcopy','off')
%         set(h,'renderer','painter');
%     else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
%     end
    print(h,[filename,'.tif'],'-dtiffn',['-r',num2str(res(7))]);
end
if FTs(8) %tiff
%     if tpb(8) %transparent background?
%         h.Color = 'none';
%         setaxbcol(AX,'none')
%         set(h,'InvertHardcopy','off')
%         set(h,'renderer','painter');
%     else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
%     end
    print(h,[filename,'.tiff'],'-dtiffn',['-r',num2str(res(8))]);
end
if FTs(9) %pdf
%     if tpb(9) %transparent background?
%         h.Color = 'none';
%         setaxbcol(AX,'none')
%         set(h,'InvertHardcopy','off')
%         set(h,'renderer','painter');
%     else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
%     end
    print(h,[filename,'.pdf'],'-dpdf');
end
if FTs(10) %png
%     if tpb(10) %transparent background?
%         h.Color = 'none';
%         setaxbcol(AX,'none')
%         set(h,'InvertHardcopy','off')
%         set(h,'renderer','painter');
%     else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
%     end
    print(h,[filename,'.png'],'-dpng',['-r',num2str(res(10))]);
end
if FTs(11) %copy
    if tpb(11) %transparent background?
        h.Color = 'none';
        setaxbcol(AX,'none')
        set(h,'InvertHardcopy','off')
        set(h,'renderer','painter');
    else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
    end
    if ispc
        print(h,'-clipboard','-dmeta') %Windows
    else
        try
            print(h,[filename,'.svg'],'-dsvg','-painters')
%             command = ['xclip -i -selection c ',filename];
            command = ['copyq write image/svg - < ',[filename,'.svg'],' && copyq select 0'];
            status = unix(command);
            delete([filename,'.svg'])
            if ~status
                error('failed.')
            end
        catch
            warning('svg-Clipboard copy fehlgeschlagen. Figure wird als pdf kopiert.')
            print(h,'-clipboard','-dpdf') %Linux / OSX
        end
    end
end
if FTs(12) %fig
    if tpb(12) %transparent background?
        h.Color = 'none';
        setaxbcol(AX,'none')
        set(h,'InvertHardcopy','off')
        set(h,'renderer','painter');
    else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
    end
    savefig([filename,'.fig'])
end
if FTs(13) %svg
%     if tpb(13) %transparent background?
%         h.Color = 'none';
%         setaxbcol(AX,'none')
%         set(h,'InvertHardcopy','off')
%         set(h,'renderer','painter');
%     else %non-transparent background?
        h.Color = bckc;
        setaxbcol(AX,col);
        set(h,'InvertHardcopy','on')
        set(h,'renderer',renderer);
%     end
    print(h,[filename,'.svg'],'-dsvg','-painters')
end
    
h.Color = bckcO;
setaxbcol(AX,colO);
set(h,'renderer',renderer);

%% Redock docked figures
h.WindowStyle = WindowStyle;

%% end of main function
end

%subfunction StrfindInCell
function [tf, idx] = StrfindInCell(C,S,option)
%StrfindInCell: Findet Strings in einem Cell-Array
%
%Syntax: [tf, idx] = StrfindInCell(C,S,option);
%              idx = StrfindInCell(C,S,option);
%   z.B. [tf, idx] = StrfindInCell(C,'This is a string',2);
%
%Eingabeparameter:
%
%   - C: Cell-Array, in dem gesucht werden soll
%   - S: Der String, nach dem gesucht werden soll
%   - option: 1 oder 2 (Double)
%           option = 1:
%           Suche nach Strings, die ganz oder teilweise im Cell vorhanden
%           sind. [z. B. gilt der String 'This is a string' als gefunden,
%           wenn nur nach dem Begriff 'string' gesucht wird]
%           option = 2:
%           Suche nur nach Cells, die komplett dem String entsprechen.
%           [z. B. wird der String 'This is a string' nicht gefunden, wenn
%           nur nach 'string' gesucht wird]
%
%Ausgabeparameter:
%
%   - idx: Index/Indizes des Cells, in denen sich der String befindet
%   - tf: True/False (logical) mit true f�r die Indizes, in denen der
%         String ganz oder teilweise vorhanden ist.
%
%Autor: Marc Jakobi 07.07.2015

if option == 1
    tf = cellfun('length',regexp(C,S)) == 1;
elseif option == 2
    tf = cellfun(@(s) ~isempty(strfind(S, s)), C);
end

if nargout == 2
    idx = find(tf);
else
    tf = find(tf); %tf gets output as idx
end
    
end


%% function to set axes background colour
function setaxbcol(AX,col)
%AX: axes object array
%col: colour
for i = 1:length(AX)
    if isstruct(col)
        AX(i).Color = col(i).col;
    else
        AX(i).Color = col;
    end
end
end

%% function for getting axes background colour
function col = getaxbcol(AX)
%col: struct containing axes background colours
col = repmat(struct,length(AX),1);
for i = 1:length(AX)
    col(i).col = AX(i).Color;
end
end
