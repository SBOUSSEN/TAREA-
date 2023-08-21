


        

clear all



np = 15;

filename='PATIENT6.mat';
filesave='PATIENT6RESULTS'

 load(filename);
    labels = fieldnames(Data);
    
     
    % Time treatement Time should be in YYYYMMDDHHmmss
    
  %Example 20230207023255

        
 
        t1 = uint64(Data.Date);
        s = num2str(t1);
         annee = str2num(s(:,1:4));
        % Extraction du vecteurs Mois
        mois = str2num(s(:,5:6));
        % Extraction du vecteur Jour
        jour = str2num(s(:,7:8));
        % Extraction du vecteur Heure
        heure = str2num(s(:,9:10));
         % Extraction du vecteur minutes
        minutes = str2num(s(:,11:12));
      
        % Extraction du vecteurs seconde
        seconde = str2num(s(:,13:14));
  

     t = datetime(annee,mois,jour,heure,minutes,seconde);
    [t,I]=sort(t,'ascend');

    

            
       
        
t=t-t(1);
Tps = day(t);

disp('Time : Ok')


     % headings name
 in = 1;
 Entete(1,in)  = {'Time'};
 
  in = in+1;
   Entete(1,in)  = {'Age'};
   in = in+1;
  Entete(1,in)  = {'GCS'};
  
 
 in = in+1;
  Entete(1,in) = {'npeaks_HR'};
 in = in+1;      
 Entete(1,in) = {'HR_HD'};
 in = in+1;
Entete(1,in) = {'HR_EHF'};    
 in = in+1;
Entete(1,in) = {'SBP_70_90'};
  in = in+1;
Entete(1,in) = {'DBP_0_50'};
in = in+1;
Entete(1,in) = {'MBP_0_65'};
in = in+1;
Entete(1,in) = {'SAPS2'};
        
  % DATA ALLOCATION
  
    HR=Data.FC;
 
    data(1,:) = unique(Tps);

    SBP=Data.PAs;
    MBP=Data.PAm;
    DBP=Data.PAd;
    GCS=Data.GCS;
    Age=Data.Age;
    SAPS2=Data.SAPS2;
    
    %GCS Transformation
    
    result=0;
    result(GCS >= 14 & GCS <= 15) = 0;
    result(GCS >= 11 & GCS <= 13) = 1;
    result(GCS >= 9  & GCS <= 10) = 2;
    result(GCS >= 6  & GCS <= 8)  = 3;
    result(GCS <= 5)              = 4;
   data(2,:)=Age;
   data(3,:)=result;
 
    
   %% HR COMPUTATION
        
           
   %HR pretreatment Mean, High Frequency and Low Frequency Content
  [HR_m,HR_v,HR_BF, HR_HF] = Moyenne(HR,np);
  % HR cutting on 24 h trunk
  sTemp = decoupage(HR_m,Tps);
  sTempHF = decoupage(HR_HF,Tps);
  sTempBF = decoupage(HR_BF,Tps);
  sTempv = decoupage(sffilt(@mean,HR_v,5),Tps);
  

 
 % HR NPeaks computation

            [a,b] = size(sTemp);
            seuilfc{1,:} = {};
            for k = 1:b
                seuilfc{1,k} = 0.25;
            end
            data(4,:) =  cellfun(@npeaks,sTempv,seuilfc);
  % Higushi computation
   HR_HFD = cellfun(@HFD_LCALC,sTemp);    
   data(5,:) = HR_HFD;
   
   %HR High frequency Conteent
            
  [HR_EHF, HR_EBF]=cellfun(@energie,sTempHF,sTempBF);
  data(6,:) = HR_EHF;
            
          
    
   %%  SBP COMPUTATION
     
           
%SBP pretreatment Mean, High Frequency and Low Frequency Content            

   [SBP_m,SBP_v,SBP_BF,SBP_HF] = Moyenne(SBP,np);
            
 % SBP cutting on 24 h trunk           
sTemp = decoupage(SBP_m,Tps);
            
  %SBP_70-90 computation     
 [a,b] = size(sTemp)

    for k = 1:b
       data(7,k) = Intervalle_REVIEW(sTemp{1,k}(1:end,1),'PAS');
             
    end
          
        disp('SBP : Ok')
        
  
        %% DBP COMPUTATION
       
%DBP pretreatment Mean, High Frequency and Low Frequency Content                      
 
            
             [DBP_m,PAD_v,DBP_BF,DBP_HF] = Moyenne(DBP,np);
          
  % DBP cutting on 24 h trunk  
         
            
            sTemp = decoupage(DBP_m,Tps);
   %DBP_0-50 computation             
            
  [a,b] = size(sTemp)          
for k=1:b
         data(8,k) = Intervalle_REVIEW(sTemp{1,k}(1:end,1),'PAD');
            
end
              
        disp('DBP : Ok')
        
        
      
    %% MBP Computation

%MBP pretreatment Mean, High Frequency and Low Frequency Content                      
             
            
             [MBP_m,MBP_v,MBP_BF,MBP_HF] = Moyenne(MBP,np);
    % MBP cutting on 24 h trunk           
            
            sTemp = decoupage(MBP_m,Tps);
            
            
 %MBP_0_65 computation 
 for k=1:b        
               data(9,k) = Intervalle_REVIEW(sTemp{1,k}(1:end,1),'PAM');
          
 end
             
   data(10,:)=SAPS2;    
           
       
        
    
    %% Conversion Mat2Cell and 
        data = num2cell(data');
        [u,v] = size(data);
        for k = 1:u
            for l = 1:v
                dStruc{1,l} = Entete{1,l};
                dStruc{k+1,l} = data{k,l};
            end
        end
        disp('Struc : Ok')  
          
        %% Results saving
        
        fName = strcat(filesave,'.mat');
     save(fName,'dStruc')
     CSV=dStruc(2:end,:);
      CSV=cell2mat(CSV);
      fName2 = strcat(filesave,'.csv');
      csvwrite(fName2,CSV);
 
 clear all
 
 %N_Peaks computation                     
       

function [npics] =  npeaks(signal,seuil)
    if (length(signal)<60)
        npics = 0;
    else
        [npics,lct] = findpeaks(signal,'MinPeakHeight',seuil);
        npics = length(npics);
    end
end

%%%%%%%%Compute Energy of the Signal

function [ EHF EBF ] = energie( XHF,YBF )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
EHF=trapz(XHF.^2);
EBF=trapz(detrend(YBF).^2);

ET=EHF+EBF;
EHF=EHF/ET;
EBF=EBF/ET;
end

%%%%%%%%%%%% Time cutting in trunk associated with the time base Temps
function [Data TT] = decoupage(signal,temps)
    Data = {};
    indice = unique(temps);
    for i = 1:length(indice)
        Temp = signal(temps == indice(i));
        Data{1,i} =Temp; %Temp(Temp~=0);
        T = temps(temps == indice(i));
        TT{1,i} = T;%(Temp~=0);
    end    
end


%%%%%%%

% Computation of the mean, variance, High et and Low frequency signal
% content
function [mF,sF,yo,y] = Moyenne(signal,nF)
signal=double(signal);
    if strcmp(inputname(1),'TEMP')
        signal(signal<20) = NaN;
    else
        signal(signal<10) = NaN;
    end
    % Suppression des pics
    mF = nanmoving_average(signal,nF,[],1); % Signal propre
    mF = smoothn(mF,'robust');
    yo = sffilt(@mean,mF,360);% Signal BF
    sF = sffilt(@var,mF,nF); % Variance
    y = mF-yo; % Signal HF
end

%%%%%%%%%%%Computation of Higushi Fractal Dimension


function [ Df,x,y,LARE ] = HFD_LCALC( data,varargin)
%HFD_LCALC calculates the Curve Lengths for various series of curves
%obtained from the input data
%   Based on Higuchi Fractal Dimension Algorithm, Df (fractal dimension)
%   can be calculated from L(k)~k^(-Df).
%   This program calculated L(k) for various k=1 to kmax
%INPUT ARGUMENTS:
%   data: Mandatory argument and should be row vector. If not, it will be
%         converted to row vector
%   kmax: If not specfied takes floor(length(data)/2) by default. If
%         specified it should be less or equal to floor(length(data)/2).
%         Else it takes the default value. (It should be an integer).


%Error check and Intialize the variables, arrays
nVar=length(varargin);
if(nVar>1)
    error('Too many Inputs::Expects only input data and optional kmax');
end
[r c]=size(data);
if (r>1)
    fprintf('Converting input data to a row vector\n');
    data=data(:)';% data needs to be a vector
end
L=length(data);
kmax=[];
if (nVar==1)
    kmax=varargin{1};
end
[r c]=size(kmax);
if((r~=1) || (c~=1) || (isnumeric(kmax)==0) || (kmax>floor(L/2)) || (kmax<1))
    kmax=floor(L/2);
    fprintf('kmax is not chosen a proper value. Changing kmax= %d\n',kmax);
end
LARE=zeros(1,kmax);%LARE stores the length of curves L(k) for k=1:kmax
y=zeros(1,kmax);%this is to store y=log(LARE)
x=zeros(1,kmax);%this is to store x=log(1/k)

for k=1:kmax
    LAk=0;
    for i=1:k
        LAi=0;
        for j=1:floor((L-i)/k)        
            LAi=LAi+abs(data(i+j*k)-data(i+(j-1)*k));
        end
        a=(L-1)/(floor((L-i)/k)*k);
        LAk=LAk+LAi*a/k;
    end
    LARE(k)=LAk/k;
    y(k)=log(LARE(k));
    x(k)=log(1/k);
    %%%% Important:: Uncomment the step i.e disp(k) while running for huge
    %%%% loops. As ctrl+c can work good to break the loop if uncommented .
    %%%% If it is commented, there is hard chance to break. Uncommenting
    %%%% will produce lot of output k=1:kmax.
    %disp(k)
end

%Df=(kmax*sum(x.*y)-sum(x)*sum(y))/(kmax*sum(x.*x)-sum(x)*sum(x));


%using polyfit 
coef=polyfit(x,y,1);
Df=coef(1);
end

%Computation of the time spent between intervalles: beware PAS=Systolic
%Blood Pressure PAM : Mean Blood Pressuree and PAD Diastolic Blood Pressure

function [plage] = Intervalle_REVIEW(signal,name)

T = 720;


    
    
    %% Pourcentage de PAs ou PAm par intervalle de temps
    if strcmp(name,'PAS')==1
        
        pl1 = find(70<signal & signal<90);
        
        if isempty(pl1) == 1
           plage = 0;
           
     
    else
        plage=100*length(pl1)/length(signal);
        end
        
        


     
    elseif strcmp(name,'PAM')==1
       
        pl1 = find(0<signal & signal<65);
        
        if isempty(pl1) == 1
            plage = 0;
       
        else
        plage=100*length(pl1)/length(signal);
        end
        
   
    elseif strcmp(name,'PAD')==1
        
        pl1 = find(0<signal & signal<50);
        
        if isempty(pl1) == 1
            plage = 0;
             else
        plage=100*length(pl1)/length(signal);
        end
        end
        
        
end
     
    
