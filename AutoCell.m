function [ROIprecomp ROIpar A B C] = AutoCell(x,ROIprecomp,ROIpar)
% function ROIprecomp = AutoCell(x)
% function [ROIprecomp ROIpar A B C] = AutoCell(x,ROIprecomp[,ROIpar])

% Input
doprepare = (nargin<2 || isempty(ROIprecomp));
douser = (nargout>1);
dodisplay = douser;
if nargin<3, ROIpar = struct; end
if ~isfield(ROIpar,'clustering'), ROIpar.clustering = []; end
if ~isfield(ROIpar,'shapepar'), ROIpar.shapepar = []; end
if ~isfield(ROIpar,'okneuron'), ROIpar.okneuron = []; end

% Data sizes
[nx ny nframes] = size(x); 
np = nx*ny;

% Process data
% (filter)
% Note that according to the way x is built in DataRegistrationScript.m,
% 1 frame represents approximately 5s (max 10s)
% So high-pass with parameter 60 corresponds approximately to 5min
% (before was 20, corresponding to 1min40s)
fn_progress('temporal high-pass','%');
xf = double(x);
for i=1:nx
    fn_progress(i/nx)
    xf(i,:,:) = fn_filt(xf(i,:,:),60,'hk',3);
end
% (remove global signal)
disp 'removing global signal'
xf = fn_imvect(xf); % 2D array space x time
xff = fn_subtract(xf, nmean(xf,2)); % subtract average frame, i.e. each pixel has average 0 value
mx = nmean(xff,1); % average signal
globcontribution = xf/mx;
xf = xf - globcontribution*mx;
xf = fn_imvect(xf,[nx ny]); % back to 3D

% STD, signals
xstd = std(xf,1,3);
csig = reshape(xf,[np nframes])';               % signals of clusters

% Create info structure to easily transfer basic variables to child
% functions
info = fn_pointer('csig',csig,'nx',nx,'ny',ny,'np',np,'rok',[]);

% Init display
if dodisplay
    initDisplay(info,x,xf,xstd);
end

% Cluster pixels
if doprepare
    ROIprecomp = precompute(info,dodisplay);
else
    %if dodisplay, set(info.hl,'ydata',ROIprecomp.cdist), end
end

% Finished the 'preparation' phase (i.e. the clustering, summarized in
% the variable tree)
if ~douser, return, end

% Interaction with user
[ROIpar clusters ineuropil alpha beta] = supervisedAutocell(info,ROIprecomp,ROIpar);

% Save preferred figure positions
prefpos = fn_get(info.figs,'position');
fn_userconfig('Autocell',prefpos)

% Finish: matrix to extract the data
% sig_neurons = (A - B*C)*sig   [here sig is np*nt]
neurons = find(ROIpar.okneuron);
nneuron = sum(ROIpar.okneuron);
nmanual = length(ROIpar.manualneurons);
% (neuropil pixels)
if isempty(ROIpar.manualneuropil)
    C = pinv(double(clusters(:)==ineuropil));
else
    C = zeros(1,np); C(ROIpar.manualneuropil) = 1/length(ROIpar.manualneuropil);
end
% (data extraction matrix)
A = sparse(nneuron+nmanual,np);
B = zeros(nneuron+nmanual,1);
for i=1:nneuron
    idxi = find(clusters(:)==neurons(i));
    alphai1 = pinv(alpha(idxi)); % inversion of sig_region -> sig
    alphai1 = alphai1 / sum(alphai1); % avoid divergences on the norm of the estimated sig_region
    A(i,idxi) = alphai1; %#ok<SPRIX>
    B(i) = sum(alphai1(:).*beta(idxi));
end
% (cancel the neuropil correction, was not such a good way to do it)
B(:)=0;
% (manual neurons)
for i=1:nmanual
    idxi = ROIpar.manualneurons{i};
    A(nneuron+i,idxi) = 1/length(idxi); %#ok<SPRIX>
end

% Close figures
close(info.figs)

%---
function initDisplay(info,x,xf,xstd)

% input
[np nx ny] = deal(info.np,info.nx,info.ny);

disp 'init display'

% define colors
col = randomcolors(np);

% display
clusters = reshape(1:np,[nx ny]);
im = reshape(col(clusters,:),[nx ny 1 3]);
% fn_figure('CLUSTERING DISTANCE')
% hl = plot(fliplr(1:np),zeros(np,1)); set(gca,'xscale','log')
% xlabel 'number of clusters'
% ylabel 'minimal distance between clusters'
% drawnow
figs = [fn_figure('STD') fn_figure('TIME COURSES') fn_figure('DATA') fn_figure('CLUSTERS')];
G = geometry; % make sure we will not link to some previous explor object present in memory
fourd(xf,G,'plot','proj',3,'in',figs(2));
a_x=fourd(x,G,'2d','shapemode','ellipse','in',figs(3));
a_xstd=fourd(xstd,G,'2d','shapemode','ellipse','in',figs(1));
a_im=fourd(im,G,'2dcol','shapemode','ellipse','in',figs(4));
set([a_x.D.cross a_xstd.D.cross a_im.D.cross],'hittest','off')
set([a_x.D a_xstd.D a_im.D],'logscale',true,'autoclipmode','prc0-0.5')

% set preferred positions
prefpos = fn_userconfig('AutoCell');
if ~isempty(prefpos), try fn_set(figs,'position',prefpos), end, end %#ok<TRYNC>

% store info
info.figs = figs;
info.a_x = a_x;
info.a_xstd = a_xstd;
info.a_im = a_im;
info.col = col;

% set user double-click -> show peak for clicked pixel
fun = @(D)userrmregion(info,[],[]);
set([a_x.D a_xstd.D a_im.D],'usercallback',fun)

% option to load collection of points to be marked (typically the centers
% of neurons from a previous experiment)
upoints = uicontrol('parent',figs(4),'string','Load points','callback',@(u,e)LoadPoints(info));
fn_controlpositions(upoints,figs(4),[1 1 0 0],[-150 -40 140 30])

    


%---
function ROIprecomp = precompute(info,dodisplay)

% input
[csig,nx,ny,np] = deal(info.csig,info.nx,info.ny,info.np);

% build connections between neighboring pixels
fn_progress('connecting pixel',np)
[ii jj] = ndgrid(1:nx,1:ny);
edges = zeros(2,2*np);
kedg = 0;
for k2=1:np
    fn_progress(k2)
    % initially every pixel is connected to its 4 neighbors
    i = ii(k2); j = jj(k2);
    if i>1
        k1 = k2-1;
        kedg = kedg+1;
        edges(:,kedg) = [k1 k2];
    end
    if j>1
        k1 = k2-nx;
        kedg = kedg+1;
        edges(:,kedg) = [k1 k2];
    end
end
edges = edges(:,1:kedg);

% replaceing NaNs with zeros
csig(isnan(csig)) = 0;

% build the graph of clusters, initialized with individual pixels
disp 'creating graph'

g = graph(csig,edges,[],[],[],'correlation');
nc = g.nv;
tree = [1:np; zeros(1,np)];
cdist = zeros(np,1);

% clustering loop
ntarget = 1;
DSTEP = 2000;
while nc>ntarget
    
    % decimate
    [tree cdist] = decimate(g,max(g.nv-DSTEP,ntarget),tree,cdist);
    nc = g.nv;
    
    % display
    fprintf('number of clusters: %i\n',nc)
    if dodisplay
        clusters = reshape(1:np,[nx ny]);
        for i=1:np, clusters(i) = clusters(tree(1,i)); end
        im = reshape(info.col(clusters,:),[nx ny 1 3]);
        info.a_im.SI.data = im;
        %     cx = reshape(g.vvalues{1}(:,clusters)',[nx ny nframes]);
        %     a_cx.SI.data = cx;
        %set(info.hl,'ydata',cdist)
        drawnow
    end
end

% delete graph, keep only the clusters information
clear g

% output
ROIprecomp = struct('tree',tree,'cdist',cdist);

%---
function [ROIpar clusters ineuropil alpha beta] = supervisedAutocell(info,ROIprecomp,ROIpar)

finished = false;
while ~finished
    % Cluster according to number of clusters selected by user
    [ROIpar.clustering clusters] = promptClustPar(info,ROIprecomp,ROIpar.clustering);
    nc = max(clusters(:));
    
    % Add manual neuron as additional clusters
    if isfield(ROIpar,'manualneurons') && strcmp(questdlg('There are manual neurons. Add them?','AutoCell','Yes','No','Yes'),'Yes')
        nmanual = length(ROIpar.manualneurons);
        for i=1:nmanual
            clusters(ROIpar.manualneurons{i}) = nc+i;
        end
        if length(ROIpar.okneuron)==nc
            ROIpar.okneuron = [ROIpar.okneuron; true(nmanual,1)];
        end
        nc = nc+nmanual;
    end
        
    % Weight pixels within regions
    %[ineuropil info.okreg alpha beta info.alpha_im info.beta_im] = initAlphaMap(info,clusters); %#ok<ASGLU>
    [ineuropil info.okreg alpha beta] = initAlphaMap(info,clusters); %#ok<ASGLU>
    
    % Compute statistics
    stats = regionStatistics(info,clusters,alpha);

    % Set user double-click -> add/remove region
    fun = @(D)userrmregion(info,clusters,alpha,stats.centers);
    set([info.a_x.D info.a_xstd.D info.a_im.D],'usercallback',fun)
    
    % Use statistics to guess which regions are actually neurons
    if isempty(ROIpar.okneuron) || length(ROIpar.okneuron)~=nc
        if isfield(ROIpar,'preset'), ROIpar.shapepar.preset = ROIpar.preset; end
        ROIpar.shapepar = promptShapePar(info,clusters,alpha,stats,ROIpar.shapepar);
        % note that promptShapePar automatically sets info.rok and shows the
        % neurons
        if isempty(ROIpar.shapepar), continue, end % window was closed, start again prompting how many clusters
    else
        info.rok = ROIpar.okneuron;
        showNeurons(info,clusters,alpha,stats.centers)
    end
    
    % Wait until user has finished modifying the selection
    hf = msgbox('Double-click on regions to select/deselect, trace new regions with right button, press OK when finished. Close figure to try a different number of clusters.','user validation');
    hu = findobj(hf,'type','uicontrol');
    set(hu,'callback',@(u,e)delete(hu));
    waitfor(hu)
    
    % Output
    finished = ishandle(hf);
    if finished
        close(hf)
        ROIpar.okneuron = info.rok;
        ROIpar.okneuron(ineuropil) = false;
    else
        ROIpar.okneuron = [];
    end
end

% Manual selections
manual = info.a_im.SI.selectionmarks;
[ROIpar.manualneuropil ROIpar.manualneurons] = deal([]);
if ~isempty(manual)
    answer = questdlg('You selected regions manually. Do you want to keep them as additional neurons?', ...
        'manual selections','Yes','No','Yes'); % ,'Yes, and use first to redefine neuropil'
    switch answer
        case 'Yes'
            ROIpar.manualneurons = {manual.dataind};
        case 'Yes, and use first to redefine neuropil'
            ROIpar.manualneuropil = manual(1).dataind;
            ROIpar.manualneurons = {manual(2:end).dataind};
    end
end

%---
function [clustpar clusters] = promptClustPar(info,ROIprecomp,clustpar)
% Decide number of clusters

% default clustering parameters
persistent defnpilclust
if isempty(defnpilclust)
    if ~isempty(clustpar) && ~isempty(clustpar.npilclust)
        defnpilclust = clustpar.npilclust; 
    else
        defnpilclust = 1000;
    end
end
if isempty(clustpar)
    clustpar = struct('nc',15000,'npilclust',[]);
end

% input
[nx ny np a_im col] = deal(info.nx,info.ny,info.np,info.a_im,info.col);

% this function can have been called after further steps have been already
% performed, then user decided to restart from the beginning, in this case
% some displays/variable need to be reset
set(info.figs(4),'name','CLUSTERS')
info.rok = [];

% loop until user confirms choice
prevpar = [];
while ~isequal(clustpar,prevpar)
    % update clusters
    clusters = reshape(1:np,[nx ny]);
    tree = ROIprecomp.tree;
    for i=find(tree(2,:)>clustpar.nc)
        clusters(i) = clusters(tree(1,i));
    end
    [u u2c c2u] = unique(clusters(:)); %#ok<ASGLU>
    clusters = reshape(c2u,[nx ny]);
    nc = clustpar.nc;
    clear u u2c c2u
    
    % cluster neuropil regions
    csiz = double(sum(sparse(clusters(:),1:nx*ny,1,nc,nx*ny),2));
    if clustpar.npilclust
        inpil = find(csiz>clustpar.npilclust);
        clusters(ismember(clusters,inpil)) = inpil(1);
        [u u2c c2u] = unique(clusters(:)); %#ok<ASGLU>
        clusters = reshape(c2u,[nx ny]);
        csiz = csiz(u);
        nc = length(u); % became probably < clustpar.nc
        clear u u2c c2u
    end
    
    % update display
    col1 = col;
    %col1(csiz<=5,:) = .2;
    a_im.SI.data = reshape(col1(clusters,:),[nx ny 1 3]);
    
    %     % put back the original signals in the time courses display
    %     cx = reshape(csig',[nx ny nframes]);
    %     a_cx.SI.data = cx;
    
    % is this choice of clustering approved?
    prevpar = clustpar;
    spec = struct('nc',{'double' 'Number of clusters'}, ...
        'npilclust','hide'); %{['xdouble [' num2str(defnpilclust) ']'] 'cluster neuropil regions larger than'});
    clustpar = fn_structedit(clustpar,spec);
    if clustpar.npilclust, defnpilclust = clustpar.npilclust; end
end
drawnow

%---
%function [ineuropil okreg alpha beta alpha_im beta_im] = initAlphaMap(info,clusters)
function [ineuropil okreg alpha beta] = initAlphaMap(info,clusters)
% (model: in every region, sig = alpha*sig_region + beta*sig_neuropil + noise)

disp 'init alpha map'

% input
[csig nx ny] = deal(info.csig,info.nx,info.ny);
nc = max(clusters(:));

% neuropil region
beta = zeros(nx,ny);
[dum ineuropil] = max(sum(sparse(clusters(:),1:nx*ny,1,nc,nx*ny),2)); %#ok<ASGLU> % a nice formula isn't it!?
idxi = find(clusters(:)==ineuropil);
sig = csig(:,idxi);
sig_neuropil = nmean(sig,2);
S = sig_neuropil;
beta(idxi) = pinv(S)*sig;
% other regions
% (first approximation: region signals are uniform average over region)
alpha = ones(nx,ny);
alpha(clusters==ineuropil) = 0;
nperclust = zeros(1,nc);
for i=1:nx*ny, ic=clusters(i); nperclust(ic)=nperclust(ic)+1; end
okreg = (nperclust>5);
okreg(ineuropil) = false;

% % prepare display
% fn_figure('WEIGHTS')
% subplot(121), alpha_im = fourd(alpha,'clip',[0 2]);
% subplot(122), beta_im  = fourd(beta);
% drawnow

% convergence loop for alpha
disp '(computation of alpha map: adjustment using a global neuropil signal might not be appropriate!)'
nloop = 3;
for kloop = 1:nloop
    fn_progress(['alpha map loop ' num2str(kloop) '/' num2str(nloop) ', region'],nc)
    for i=1:nc
        if ~okreg(i), continue, end % note that ineuropil will not go through
        fn_progress(i)
        idxi = find(clusters(:)==i);
        sig = csig(:,idxi);
        alphai1 = pinv(alpha(idxi)); % inversion of sig_region -> sig
        alphai1 = alphai1 / sum(alphai1); % avoid divergences on the norm of the estimated sig_region
        sig_region = (sig - sig_neuropil*beta(idxi)')*alphai1';
        S = [sig_region sig_neuropil];
        ab = S\sig;
        % no negative contribution of neuropil authorized
        negnpil = ab(2,:)<0;
        ab(2,negnpil) = 0;
        S = sig_region;
        ab(1,negnpil) = S\sig(:,negnpil);
        alpha(idxi) = ab(1,:);
        beta(idxi) = ab(2,:);
    end
    %     % display
    %     alpha_im.SI.data = alpha;
    %     beta_im.SI.data  = beta;
end

%---
function stats = regionStatistics(info,clusters,alpha)

% Compute regions statistis
disp 'regions statistics'
[ii jj] = ndgrid(1:info.nx,1:info.ny);
% size and 'roundness' of regions
nc = max(clusters(:));
[csiz dispersion csiz2 border] = deal(zeros(nc,1));
centers = zeros(2,nc);
for i=find(info.okreg)
    % pixels inside the region
    maski = (clusters==i);
    idxi = find(maski);
    % square root of number of pixels, represents "length" of the region
    csiz(i) = sqrt(length(idxi));
    % "dispersion" from the center, gives small number for round shapes,
    % high number for elongated shape
    ic = mean(ii(idxi));
    jc = mean(jj(idxi));
    devi = sqrt(mean((ii(idxi)-ic).^2) + mean((jj(idxi)-jc).^2));
    dispersion(i) = devi/csiz(i);
    % weighted center and dispersion
    ai = alpha(idxi);
    a = sum(ai);
    ic = sum(ii(idxi).*ai)/a;
    jc = sum(jj(idxi).*ai)/a;
    csiz2(i) = sqrt(sum((ii(idxi)-ic).^2.*ai)/a + sum((jj(idxi)-jc).^2.*ai)/a);
    % "border": ratio of average value on the region border to average
    % value on the full region
    
    % replace 'borderi = bwmorph(maski,'remove');' by the lines below to
    % avoid using the Image Toolbox
    borderi = maski & ~(maski(:,[2:end end]) & maski(:,[1 1:end-1]) & maski([2:end end],:) & maski([1 1:end-1],:));
    
    valglob = mean(alpha(maski));
    valborder = mean(alpha(borderi));
    border(i) = valborder/valglob;
    % center of cell
    centers(:,i) = [ic; jc];
end
stats = struct('csiz',csiz,'dispersion',dispersion,'csiz2',csiz2,'border',border, ...
    'centers',centers);

%---
function showNeurons(info,clusters,alpha,centers)

% color selected neurons brighter
col1 = info.col;
bad = ~info.rok;
col1(~info.rok,:) = info.col(bad,:)/5;
col1(~info.okreg,:) = .2; % neuropil and <=5 pixels regions
colclusters = reshape(col1(clusters,:),[info.nx info.ny 1 3]);
colclusters = fn_mult(colclusters,alpha);
info.a_im.SI.data = colclusters;

% also indicate selected neurons in the other displays
pp = centers(:,info.rok)-1;
deco = selectionND('point2D',pp);
[info.a_x.SI.decoration info.a_xstd.SI.decoration] = deal(deco);


%---
function shapepar = promptShapePar(info,clusters,alpha,stats,shapepar)

% default shape criteria
if isempty(shapepar), shapepar = struct; end
if isfield(shapepar,'preset')
    if strfind(shapepar.preset, '10X obj')
        minsize = 4;
    else
        minsize = 8;
    end
    shapepar = rmfield(shapepar,'preset');
else
    minsize = 8;
end
defshapepar = struct('minsize',minsize,'maxsize',50,'maxdispersion',.5,'maxborder',.9);
shapepar = fn_structmerge(defshapepar,shapepar,'strict');

% loop until user confirms choice
prevpar = [];
while ~isequal(shapepar,prevpar)
	info.rok = (stats.csiz>=shapepar.minsize & stats.csiz<=shapepar.maxsize ...
        & stats.dispersion<=shapepar.maxdispersion & stats.border<=shapepar.maxborder);
    showNeurons(info,clusters,alpha,stats.centers)
    prevpar = shapepar;
    %disp([num2str(sum(info.rok)) ' regions selected'])
    set(info.figs(4),'name',sprintf('%i REGIONS',sum(info.rok)))
    shapepar = fn_structedit(shapepar);
    if isempty(shapepar), return, end
end
    
%---
function userrmregion(info,clusters,alpha,centers)

% input
[D csig rok] = deal(info.a_im.D,info.csig,info.rok);

% select time point of the peak
ij = D.SI.ij;
tc = csig(:,sub2ind(D.SI.sizes,ij(1),ij(2)));
[m tidx] = max(tc); %#ok<ASGLU>
talreadyset = (D.SI.G.ijkl2(3) == tidx);
if ~talreadyset
    D.SI.G.ijkl2(3) = tidx;
end

% stop here if rok is not valid
if isempty(rok), return, end

% toggle the state of the selected region
kcluster = clusters(ij(1),ij(2));
if ~info.okreg(kcluster), return, end
curstate = rok(kcluster);
if curstate && ~talreadyset, return, end
rok(kcluster) = ~curstate;
disp([num2str(sum(rok)) ' regions selected'])
set(info.figs(4),'name',sprintf('%i REGIONS',sum(rok)))
info.rok = rok;

% update display
showNeurons(info,clusters,alpha,centers)


%---
function LoadPoints(info)

% remove all currently marked points
delete(findobj(0,'tag','loadedpoints'))

% previous segmatation result
f = fn_getfile('*regions.mat','Select file with previous segmentation result');
if isequal(f,0), return, end
[A nx ny avgimg] = fn_loadvar(f,'A','nx','ny','avgimg');
name = fn_fileparts(fn_fileparts(f,'path'),'base');

% show average image and segmentation result from last experiment
fn_figure(name)
colormap(gray(256))
ha1 = subplot(121);
a4d = fourd(avgimg,'2d','clip',fn_clip(avgimg','prc.1','getrange'));
Gprev = a4d.G; % geometry for this display
% imagesc(fn_clip(avgimg','prc.1'),[0 1]), axis image
ha2 = subplot(122);
a = reshape(sum(fn_mult(A,sum(A>0,2)),1),[nx ny]);
fourd(a,'2d',Gprev,'clip',[0 2])
% imagesc(a',[0 2]), axis image

% detect cell centers (simply take the brigtest pixel of each region)
[ii jj] = ndgrid(1:nx,1:ny);
ii = column(ii); jj = column(jj);
ncell = size(A,1);
ipos = zeros(2,ncell);
At = A'; % will make 'find' below much faster
for i=1:ncell
    [idx, ~, val] = find(At(:,i));
    ipos(1,i) = sum(ii(idx).*val);
    ipos(2,i) = sum(jj(idx).*val);
end


% display
axeshandles = [info.a_x.D.ha info.a_im.D.ha info.a_xstd.D.ha];
nwithmotion = length(axeshandles);
axeshandles = [axeshandles ha1 ha2]; % rotation/translation will not be applied to points in these displays
for k=1:length(axeshandles)
    ha = axeshandles(k);
    hl(k)=line(ipos(1,:)-1,ipos(2,:)-1,'parent',ha,'tag','loadedpoints','hittest','off', ...
        'linestyle','none','color','r','marker','+','markersize',4);
end
hlprev = hl(1:nwithmotion);
hlcur = hl(nwithmotion+1:end);

% controls for translation and rotation
s = struct('x',0,'y',0,'theta',0);
fn_control(s,@(s)ApplyMotion(s,ipos,nx,ny,hlcur,hlprev,Gprev))

%---
function ApplyMotion(s,ipos,nx,ny,hlprev,hlcur,Gprev)

% build rotation/translation matrix
T1 = [1 0 0; [-nx/2; -ny/2] eye(2)]; 
R = [1 0 0; [0; 0] [cosd(s.theta) -sind(s.theta); sind(s.theta) cosd(s.theta)]]; 
T2 = [1 0 0; [s.x+nx/2-1; s.y+ny/2-1] eye(2)]; 
mat = T2*R*T1;

% apply to the geometry object
Gprev.mat = mat;
set(hlprev,'xdata',ipos(1,:)+Gprev.grid(1,2),'ydata',ipos(2,:)+Gprev.grid(2,2))

% apply to points that moved
xx = [ones(1,size(ipos,2)); ipos];
ipos = mat(2:3,:) * xx;
set(hlcur,'xdata',ipos(1,:),'ydata',ipos(2,:))

% points that did not move: coordinate system might have changed however!!!







