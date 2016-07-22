function stark_split5

warning off;
ppr_size = [14.8 11.8];

% Initialisierung

e0 = 8.854187817e-12;

g = 3e-3;
r = 15e-3/2;

d_glas = 0.7e-3;
d_al2o3 = 0.2e-3;
d_bso = 0.7e-3;

e_glas = 7.6;
e_al2o3 = 10.55;
e_bso = 56;

A = pi*(r)^2;

c_ext = 1e-9;
c_gap = e0*A/g;
c_glas = e_glas*e0*A/d_glas;
c_al2o3 = e_al2o3*e0*A/d_al2o3;
c_bso = e_bso*e0*A/d_bso;

c_diel = 1/(1/c_bso+1/c_glas+1/c_al2o3);

c_bd = (c_gap*c_diel)/(c_gap+c_diel);


% Dateipfade der Messdaten.

loc_dat = 'D:\documents\git\ag_praktikum_2016\2016-07-01\Daten\';
loc_main = 'D:\documents\git\ag_praktikum_2016\2016-07-01\';

% loc_dat= '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-22/Daten';
% loc_main = '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-22';

% Beschaffe mir die Dimensionen.

cd(loc_dat);
tmp = importdata('16Jul00001_RTO.dat');
tmp = tmp.data;
[n channel] = size(tmp);
time_volt = tmp(:,1);
time_delta = time_volt(2)-time_volt(1);
cd(loc_main);

% Datei-Nummern und Wellenlängen.

m=1620;
data_nmb = 1:1:m;
loop_nmb = 1:1:15;
wvl_nmb = 1:1:36;
vert_nmb = 1:1:3;

raw_data = zeros(n,m);

rogowski = zeros(n,m);

wavelength = 491.8:0.02:492.5;
vertical = 6.8:0.15:7.1;

% Entnehmen der Datein.
% Es kommen zuerst die 10 loops zu einer Wellenlänge. Dann erneut 10 loops
% für die nächste. Macht 360 Datein, bis eine Höhe abgearbeitet ist.

cd(loc_dat);

for h=vert_nmb-1;
    
    for w=wvl_nmb-1;
        
        for l=loop_nmb;

            file_nmb = num2str(l+15*w+540*h,'%05d');
            nmb = l+15*w+540*h;
            disp(file_nmb);
            file = strcat('16Jul',file_nmb,'_RTO.dat');
            tmp = importdata(file);
            tmp = tmp.data;

            raw_data(:,nmb) = tmp(:,channel-1);
            
            rogowski(:,nmb) = tmp(:,channel);
        
        end
        
    end
    
end

cd(loc_main);
save('base_dat.mat');
% load('base_dat.mat');

% Kontrolle des Experimentes
    
    f = figure;
    hold on;
    meshc(time_volt/1e-6,linspace(1,1620,1620),rogowski');view(2);
    xlabel('time in µs');
    ylabel('vertical slit pos. in inch');
    c = colorbar;
    c.Label.String = 'voltage in V';
    title('rowgowski coil via experiment duration');
    set(gca,'YTick',[270 810 1350]);
    set(gca,'YTickLabel',{'6.8' '6.95' '7.1'});
    
    set(c,'fontsize',12);
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'rowgowski_full','pdf');
%     print('rowgowski_full2','-dpdf','-noui','-bestfit');
    
%     savefig('rowgowski_full');
    hold off; close(f);
    
% FFT abziehen der ersten 200 Punkte je Messung.

for i=1:1620  
    
    tmp = raw_data(1:270,i);
    tmpvec = [tmp' tmp' tmp' tmp' tmp' tmp'];
    tmpvec = tmpvec(1:n);
    tmpfft = fft(tmpvec);
    raw_fft(:,i) = fft(raw_data(:,i))-tmpfft';
    raw_fft(10:end-10,i) = 0;
    raw_data(:,i) = real(ifft(raw_fft(:,i)));
    
end

for i=1:36;
    
    stark_68in(:,i) = 1/15*sum(raw_data(:,(i-1)*15+1:(i-1)*15+15),2);
    stark_695in(:,i) = 1/15*sum(raw_data(:,(i-1)*15+540:(i-1)*15+555),2);
    stark_71in(:,i) = 1/15*sum(raw_data(:,(i-1)*15+1080:(i-1)*15+1095),2);
    
end

%% Bilder der Stark-Aufspaltungen.

% Wellenlängen-Offset. SGOLAY-Filt.
% Abziehen von Wellenlängen, bei denen keine Emission stattfindet.

    wvl68 = 1/8*sum(stark_68in(:,1:8),2);
    wvl695 = 1/8*sum(stark_695in(:,1:8),2);
    wvl71 = 1/8*sum(stark_71in(:,1:8),2);

    for i=1:36
        stark68_korr(:,i) = stark_68in(:,i)-wvl68;
        stark695_korr(:,i) = stark_695in(:,i)-wvl695;
        stark71_korr(:,i) = stark_71in(:,i)-wvl71;
    end
    
% SGOLAY-Filt.

    sgf_71in = sgolayfilt(stark71_korr,2,51);
    sgf_695in = sgolayfilt(stark695_korr,2,51);
    sgf_68in = sgolayfilt(stark68_korr,2,51); 
    
% Bilder der Stark-Aufspaltungen.

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,sgf_68in');view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
c = colorbar;
c.Label.String = 'intensity, a.u.'; 
title('vertical position 6.8 inch');

    set(c,'fontsize',12);
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_68in','pdf');
%     print('stark_68in2','-dpdf','-noui','-bestfit');
 
% savefig('stark_68in.fig');
hold off; close(f);

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,sgf_695in');view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
title('vertical position 6.95 inch');c = colorbar;
c.Label.String = 'intensity, a.u.';    

    set(c,'fontsize',12);
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_695in','pdf');
%     print('stark_695in2','-dpdf','-noui','-bestfit');

% savefig('stark_695in.fig');
hold off; close(f);

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,sgf_71in');view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
title('vertical position 7.1 inch');c = colorbar;
c.Label.String = 'intensity, a.u.';   

    set(c,'fontsize',12);
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_71in','pdf');
%     print('stark_71in2','-dpdf','-noui','-bestfit');
    
% savefig('stark_71in.fig');
hold off; close(f);

% Schaue mir die Entladungscharakteristik an.

for i = [6.8 6.95 7.1]
    
    cd(loc_dat);
    h = num2str(i);
    
    file = strcat('16Jun_electrics_',h,'mm_RTO.dat');
    tmp = importdata(file);
    cd(loc_main);

    tmp = tmp.data;
    chrg_tmp(:) = tmp(:,3);
    volt_appl(:) = tmp(:,2);
    volt_diff(:) = 1/time_delta*diff(tmp(:,2));
    chrg_diff(:) = 1/time_delta*diff(tmp(:,3));
   
chrg_tmp = chrg_tmp*c_ext;
chrg_diff = c_ext*chrg_diff;
    
% Fit des C_tot und Berechnung des Entladungsstromes.

    xdata = volt_appl(1:150);
    ydata = chrg_tmp(1:150);
    
    P = polyfit(xdata,ydata,1);
    c_tot = P(1);
    
    c_par = c_tot-c_bd;
    
    volt_gap = volt_appl*(1-c_par/c_diel)-chrg_tmp./c_diel;
    current_dis = (1+c_gap/c_diel)*(chrg_diff-c_tot*volt_diff);

% Check.

    lname1 = strcat('lissajous',h,'.bmp');
    lname2 = strcat('lissajous2',h,'.bmp');
    
    cname1 = strcat('currentdis',h,'.bmp');
    cname2 = strcat('currentdis2',h,'.bmp');

    f = figure;hold on;
    plot(volt_appl,chrg_tmp);
    xlabel('U_{appl}/V');
    ylabel('Q_{ext}/C');
    title('applied voltage over total charge');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,lname1,'bmp');
%     print(lname2,'-dpdf','-noui','-bestfit');
%     savefig(sprintf('lissajous%s.fig',h));
    hold off; close(f);
    
    f = figure;
    hold on;
    plot(time_volt(1:999)/1e-6,volt_gap(1:999),'r',time_volt(1:999)/1e-6,volt_appl(1:999),'r-.');
    xlabel('time in µs');
    yyaxis left
    ylabel('voltage in V');
    
    plot(time_volt(1:999)/1e-6,current_dis(1:999)*100000,'c');
    yyaxis right
    ylabel('current in mA');
    
    legend('U_{gap}','U_{app}','I_{dis}');
    title('current/appl. & gap voltage via time');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,cname1,'bmp');
%     print(cname2,'-dpdf','-noui','-bestfit');
%     savefig(sprintf('currentdis%s.fig',h));
    hold off; close(f);
    
end
    
% Mache die Stark-Aufspaltung.

clear trsh Ind I inmax max nmax

[trsh, Ind] = max(abs(stark71_korr));
inmax = 0;

    for i=1:36;
        
        j=Ind(i);
        nmax = stark71_korr(j,i);
        
        if (abs(inmax)<=abs(nmax))
            nInd = [i j];
            inmax = nmax;
        end
    end
    
    time = time_volt(nInd(2))*1e6;
    
    in1max = 0;
    in2max = 0;
    starkshift = -stark71_korr(nInd(2),:);
    
    for i=1:14;
        
        tmp1 = starkshift(1+i);
        tmp2 = starkshift(15+i);
        
        if (abs(in1max)<=max(tmp1))
            wInd1=i+1;
            in1max = tmp1;
        end
        
        if (abs(in2max)<=max(tmp2))
            wInd2=i+15;
            in2max = tmp2;
        end
        
    end
            
    sep = wavelength(wInd2)-wavelength(wInd1);
    
    fieldstrengthsq = -58.557+18.116*sep+3130.96*(sep)^2+815.6*(sep)^3;
    fieldstrength = sqrt(fieldstrengthsq);
    
    f = figure;hold on;
    plot(wavelength,starkshift);
    xlabel('wavelength in nm');
    ylabel('intensity, a.u.');
    title(sprintf('shift, 7.1in, %gµs, E=%g kV/cm', time, fieldstrength));
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_shift71in','bmp');
%     print(lname2,'-dpdf','-noui','-bestfit');
%     savefig(sprintf('lissajous%s.fig',h));
    hold off; close(f);


%% Daten.

save('2016-07-01.mat');

end