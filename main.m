% This script can be used for the decomposition of LC-DAD data, in oder to
% obtain the individual UV-vis spectrum for each substance detected in MS.
% Overlapping spectra are resolved by a cnmf (Constrained Non-negative
% Matrix Factorization) model, with the power of UHPLC-DAD-UHRMS (Ultra
% High Performance Liquid Chromatography-Diode Array Detector-Ultra High
% Resolution Mass Spectrometry).

% For more infomation please refer to our paper:
% Xu N, Hu M, Li X et al. Resolving ultraviolet-visible spectral for
% complex dissolved mixtures of multitudinous organic matters in aerosol.
% Analytical Chemistry.

%% inputting variants
clear
clc
warning off
load data peakid peakFWHM peakms peakrt
% data=data; - struct array with fields: name(optional), DAD, time, and peaks. The row number equals to the number of samples.
% data.name - char vector of the sample name
% data.DAD - numeric array of the LC-DAD data matrix of the sample. Rows represent the wavelength, while columns represent the retention time.
% data.time - cell of extracted ion chromatography time vector for each peak in the sample. The column number of the cell euquals to the peak number.
% data.peaks - cell of extracted ion chromatography intensity vector for each peak in the sample. The column number of the cell euquals to the peak number.
% peakid - numeric vector of the identification number of peaks. The column number euquals to the peak number.
% peakFWHM - numeric vector of the full width at half maximum of peaks. The column number euquals to the peak number.
% peakms - numeric vector of the m/z detected by MS of peaks. The column number euquals to the peak number.
% peakrt - numeric vector of the retention times of peaks. The column number euquals to the peak number.



%% set parameters
%time unit-min; wavelength unit-nm;
Info.tc=300; % retention time sampling points per minute
Info.te=0.09; % retention time deviation between DAD and MS
Info.tstart=0; % start time of LC-DAD data for decomposition
Info.tend=90; % end time of LC-DAD data for decomposition
Info.tinterval=0.1; % moving time step of the window
Info.twindow=2; % the time width of the window
Info.wlc=1; % sampling interval of wavelength
Info.wl1=250; % start wavelength of LC-DAD data for decomposition
Info.wl2=700; % end wavelength of LC-DAD data for decomposition

% calculating parameters
cw=round(Info.twindow*Info.tc);
Info.twindow=cw/Info.tc;
ci=round(Info.tinterval*Info.tc);
Info.tinterval=ci/Info.tc;
cs=round(Info.tstart*Info.tc);
Info.tstart=cs/Info.tc;
t1=cs/Info.tc;
ce=round(Info.tend*Info.tc);
Info.tend=ce/Info.tc;
c2=min(cs+cw,ce);
t2=c2/Info.tc;
peakrt=peakrt-Info.te;

%% initializing window
samplenum=size(data,2);
nc=zeros(1,samplenum);
a=cell(1,samplenum);
an=cell(1,samplenum);
status=false(1,samplenum);
retc=zeros(1,samplenum);
Height=zeros(1,samplenum);
FWHM=zeros(1,samplenum);
Spectrum=zeros(1,(Info.wl2-Info.wl1)*Info.wlc+1);
Spectrum(1,:)=Info.wl1:1/Info.wlc:Info.wl2;
result=cell(1,6);
result{1,1}='No.';
result{1,2}='rt';
result{1,3}='id';
result{1,4}='ms';
result{1,5}='period';
result{1,6}='hr';
id=2;

% extracting and splicing DAD data among samples with normalization
for ii=1:samplenum
    a{ii}=data(ii).DAD(:,cs+1:c2);
    nc(ii)=norm(max(a{ii},0.2));
    an{ii}=a{ii}/nc(ii);
end
aa=cell2mat(an);

%% moving window
while t2<=Info.tend
    % initialization
    [n,m]=size(a{1});
    % finding peaks of which the retention time lies in the window
    MSid=find(peakrt+max(peakFWHM,0.15)*1>t1 & peakrt-max(peakFWHM,0.15)*1<t2);
    sl = length(MSid);
    if sl~=0
        rmax=min(samplenum*m,n);
        if sl>rmax
            sl=rmax;
            MSid=MSid(1:rmax);
        end
        R=sl;
        % extracting and splicing MS data of extracted ion chromatography of each found peak among samples
        h0=cell(R,samplenum);
        h02=cell(R,samplenum);
        for j = 1:R
            jid=MSid(j);
            for ii=1:samplenum
                hs2=interp1(data(ii).time{jid}-Info.te,data(ii).peaks{jid},t1+1/Info.tc:1/Info.tc:t2,'makima');
                hs=hs2/nc(ii);
                h0{j,ii}=hs;
                h02{j,ii}=hs2;
            end
        end
        % obtaining normalized initial eluting peak matrix
        H0=cell2mat(h0);
        H02=cell2mat(h02);
        H0(H0<0)=0;
        ex=find(max(H02,[],2)-min(H02,[],2)<20000);
        H0(ex,:)=[];
        MSid(ex)=[];
        R = length(MSid);
        % calculating corresponding initial spectrographic matrix
        s=svd(H0);
        w0=aa*pinv(H0,prctile(s,20));
        w0=max(w0,0).^1;
        w0=max(smoothdata(w0,1,'lowess',n/10),0);
        w0=w0+max(w0)/10;
        
        % iteration
        % calling cnmf function to do non-negative matrix factorization
        [wR,hR] = cnmf(aa,R,samplenum,'W0',w0,'H0',H0,'Algorithm','mult','Options',statset('MaxIter',100,'TolFun',0.0001,'TolX',0.0001));%'W0',w0,,'UseParallel',true
        
        % interpretation
        % checking whether to output the peak
        for j = 1:R
            jid=MSid(j);
            if  peakrt(jid)<t2-max(peakFWHM(jid)*1,0.15) && peakrt(jid)>t1+max(peakFWHM(jid)*1,0.15)
                if  max(wR(:,j))==0
                    peakms(jid)=-9999;
                else
                    statusall=0;
                    hRR=zeros(samplenum,m);
                    for ii = 1:samplenum
                        hRR(ii,:)=hR(j,(ii-1)*m+1:ii*m)*nc(ii);
                    end
                    if max(hRR,[],'all')<0.01
                        peakms(jid)=-9999;
                    else
                        mh=mean(hRR);
                        mmh=max(mh);
                        if abs(mh(m-6)-mh(m-1))<mmh/20 && mh(m-1)<mmh/20 && mh(2)<mmh/20
                            for ii = 1:samplenum
                                hn=hRR(ii,:);
                                status(ii)=abs(hn(m-6)-hn(m-1))<max(max(hn)/100,0.005) && hn(m-1)<max(max(hn)/100,0.01) && hn(2)<max(max(hn)/100,0.01);
                                statusall=statusall+status(ii);
                            end
                            
                            if statusall>samplenum*0.95
                                % recognized as a complete eluting peak in the window
                                for ii = 1:samplenum
                                    [~,Ih]=max(hRR(ii,:));
                                    retc(ii)=t1+Ih/Info.tc;
                                    a{ii}=max(a{ii}-wR(:,j)*hRR(ii,:),0); % subtracted from the DAD data matrix
                                    [Height(ii),~]=max(hRR(ii,:));
                                    FWHM(ii)=sum(hRR(ii,:)>Height(ii)/2)/Info.tc;
                                end
                                result{id,1}=id-1;
                                result{id,2}=mean(retc);
                                % outputting eluting peak
                                result{id,5}=[t1,t2];
                                result{id,6}=hRR;
                                % outputting absorbing spectrum
                                Spectrum(id,:)=wR(:,j);
                                % outputting peak identification
                                if abs(result{id,2}-peakrt(jid))<max(peakFWHM(jid),0.15) *1.5 && mean(rmoutliers(FWHM,'gesd'))<max(peakFWHM(jid),0.15) *1.5
                                    result{id,3}=peakid(jid);
                                    result{id,4}=peakms(jid);
                                else
                                    result{id,3}='unknown';
                                    result{id,4}='unknown';
                                end
                                id=id+1;
                                peakms(jid)=-9999;
                            end
                        end
                    end
                end
            end
        end
        
        % removing peaks not participated in the following calculation
        ca=find(peakms==-9999);
        peakms(ca,:)=[];
        peakrt(ca,:)=[];
        peakFWHM(ca,:)=[];
        peakid(ca,:)=[];
        for ii=1:samplenum
            data(ii).time(ca,:)=[];
            data(ii).peaks(ca,:)=[];
        end
        
    end
    
    % checking if the window has moved to the end
    if t2==Info.tend
        break
    end
    
    % extracting and splicing DAD data among samples for next window with normalization
    for ii=1:samplenum
        a{ii}=[a{ii},data(ii).DAD(:,c2+1:min(ce,c2+ci))];
        a{ii}=a{ii}(:,end-cw+1:end);
        nc(ii)=norm(max(a{ii},0.2));
        an{ii}=a{ii}/nc(ii);
    end
    aa=cell2mat(an);
    aa(aa<=0)=0;
    c2=min(ce,c2+ci);
    t2=c2/Info.tc;
    t1=t2-Info.twindow;
end
