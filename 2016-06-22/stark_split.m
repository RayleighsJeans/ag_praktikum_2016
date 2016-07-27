function stark_split

warning off;

%% Initialisierung

% e0 = 8.854187817e-12;
% 
% g = 3e-3;
% r = 15e-3/2;
% 
% d_glas = 0.7e-3;
% d_al2o3 = 0.2e-3;
% d_bso = 0.7e-3;
% 
% e_glas = 7.6;
% e_al2o3 = 10.55;
% e_bso = 56;
% 
% A = pi*(r)^2;
% 
% c_ext = 1e-9;
% c_gap = e0*A/g;
% c_glas = e_glas*e0*A/d_glas;
% c_al2o3 = e_al2o3*e0*A/d_al2o3;
% c_bso = e_bso*e0*A/d_bso;
% 
% c_diel = 1/(1/c_bso+1/c_glas+1/c_al2o3);
% 
% c_bd = (c_gap*c_diel)/(c_gap+c_diel);


%% Dateipfade der Messdaten.

loc_dat = 'D:\documents\git\ag_praktikum_2016\2016-06-22\Daten';
loc_main = 'D:\documents\git\ag_praktikum_2016\2016-06-22\';

% loc_dat= '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-22/Daten';
% loc_main = '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-22';

%% Beschaffe mir die Dimensionen.

% cd(loc_dat);
% tmp = importdata('16Jun35001_RTO.dat');
% tmp = tmp.data;
% [n channel] = size(tmp);
% time_volt = tmp(:,1);
% time_delta = time_volt(2)-time_volt(1);
% cd(loc_main);

%% Datei-Nummern und Wellenlängen.

% m=1440;
% data_nmb = 1:1:m;
% loop_nmb = 1:1:10;
% wvl_nmb = 1:1:36;
% vert_nmb = 1:1:4;
% 
% raw_data = zeros(n,m);
% 
% rogowski = zeros(n,m);
% 
% wavelength = 491.8:0.02:492.5;
% vertical = 6.8:0.1:7.1;

%% Entnehmen der Datein.
%Es kommen zuerst die 10 loops zu einer Wellenlänge. Dann erneut 10 loops
%für die nächste. Macht 360 Datein, bis eine Höhe abgearbeitet ist.

% cd(loc_dat);
% 
% for h=vert_nmb-1;
%     
%     for w=wvl_nmb-1;
%         
%         for l=loop_nmb;
% 
%             file_nmb = num2str(35000+l+10*w+360*h);
%             nmb = l+10*w+360*h;
%             disp(file_nmb);
%             file = strcat('16Jun',file_nmb,'_RTO.dat');
%             tmp = importdata(file);
%             tmp = tmp.data;
% 
%             raw_data(:,nmb) = tmp(:,4)+0.16;
%             
%             rogowski(:,nmb) = tmp(:,5);
%         
%         end
%         
%     end
%     
% end

cd(loc_main);
% save('base_dat.mat');

% Skippe die Schleifen und Initialiserung, lade Basis-Datensatz.
load base_dat.mat

% Kontrolle des Experimentes
    
    f = figure;
    hold on;
    meshc(time_volt/1e-6,linspace(1,1440,1440),rogowski'+0.160);view(2);
    xlabel('time in µs');
    ylabel('vertical slit pos. in inch');
    c = colorbar;
    c.Label.String = 'voltage in V';
    % title('rogowski coil via experiment duration');
    set(gca,'YTick',[180 540 900 1260]);
    set(gca,'YTickLabel',{'6.8' '6.9' '7.0' 'kathode'});
    ax = gca;
    ax.YTickLabelRotation = 90;
    box on;set(gca,'Layer','top');
    
%     saveas(gcf,'rogowski_full','pdf');
    print('rogowski_full','-dpdf','-noui','-bestfit');
    
%     savefig('rogowski_full');
    hold off; close(f);
    
% FFT abziehen der ersten 200 Punkte je Messung.

for i=1:1440   
    
    tmp = raw_data(1:150,i);
    tmpvec = [tmp' tmp' tmp' tmp' tmp' tmp' tmp'];
    tmpvec = tmpvec(1:1000);
    tmpfft = fft(tmpvec);
    raw_fft(:,i) = fft(raw_data(:,i))-tmpfft';
    raw_fft(10:end,i) = 0;
    raw_data(:,i) = real(ifft(raw_fft(:,i)));
    
end
    
% Einzelne Höhen herausfilter. Mache Mittelung über 10 loops. 

for i=1:36;
    
    stark_68in(:,i) = 1/10*sum(raw_data(:,(i-1)*10+1:(i-1)*10+10),2);
    stark_69in(:,i) = 1/10*sum(raw_data(:,(i-1)*10+361:(i-1)*10+370),2);
    stark_7in(:,i) = 1/10*sum(raw_data(:,(i-1)*10+721:(i-1)*10+730),2);
    stark_71in(:,i) = 1/10*sum(raw_data(:,(i-1)*10+1081:(i-1)*10+1090),2);
    
end


%% Wellenlängen-Offset. SGOLAY-Filt.

% Abziehen von Wellenlängen, bei denen keine Emission stattfindet.

    wvl68 = 1/4*sum(stark_68in(:,1:4),2);
    wvl69 = 1/4*sum(stark_69in(:,1:4),2);
    wvl7 = 1/4*sum(stark_7in(:,1:4),2);
    wvl71 = 1/4*sum(stark_71in(:,1:4),2);

    for i=1:36
        stark68_korr(:,i) = stark_68in(:,i)-wvl68;
        stark69_korr(:,i) = stark_69in(:,i)-wvl69;
        stark7_korr(:,i) = stark_7in(:,i)-wvl7;
        stark71_korr(:,i) = stark_71in(:,i)-wvl71;
    end

    sgf_71in = sgolayfilt(stark71_korr,2,151);
    sgf_7in = sgolayfilt(stark7_korr,2,151);
    sgf_69in = sgolayfilt(stark69_korr,2,151);
    sgf_68in = sgolayfilt(stark68_korr,2,151); 
    
% Bilder der Stark-Aufspaltungen.

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,-sgf_68in'-min(min(-sgf_68in)));view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
% title('vertical position 6.8 inch');
c = colorbar;
c.Label.String = 'intensity, a.u.';
box on;set(gca,'Layer','top');

%     saveas(gcf,'stark_68in','pdf');
    print('stark_68in','-dpdf','-noui','-bestfit');
 
% savefig('stark_68in.fig');
hold off; close(f);

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,-sgf_69in'-min(min(-sgf_69in)));view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
% title('vertical position 6.9 inch');
c = colorbar;
c.Label.String = 'intensity, a.u.';
box on;set(gca,'Layer','top');

%     saveas(gcf,'stark_69in','pdf');
    print('stark_69in','-dpdf','-noui','-bestfit');
 
% savefig('stark_69in.fig');
hold off; close(f);

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,-sgf_7in'-min(min(-sgf_7in)));view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
% title('vertical position 7 inch'); 
c = colorbar;
c.Label.String = 'intensity, a.u.';
box on;set(gca,'Layer','top');

%     saveas(gcf,'stark_7in','pdf');
    print('stark_7in','-dpdf','-noui','-bestfit');
 
% savefig('stark_7in.fig');
hold off; close(f);

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,-sgf_71in'-min(min(-sgf_71in)));view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
% title('vertical position 7.1 inch');
c = colorbar;
c.Label.String = 'intensity, a.u.';
box on;set(gca,'Layer','top');
    
%     saveas(gcf,'stark_71in','pdf');
    print('stark_71in','-dpdf','-noui','-bestfit');
 
% savefig('stark_71in.fig');
hold off; close(f);

% Schaue mir die Entladungscharakteristik an.

for i = [6.8 6.9 7.0 7.1]
    
    cd(loc_dat);
    h = num2str(i);
    
    if (i==7)
        h = '7.0';
    else
        
    end
    
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

    lname = strcat('lissajous',h,'.bmp');
    
    cname = strcat('currentdis',h,'.bmp');

    f = figure;hold on;
    plot(volt_appl,chrg_tmp);
    x = volt_appl;
    y = chrg_tmp;    
    axis([min(x) max(x) min(y)-0.1*max(abs(y)) max(y)+0.1*max(abs(y))]);
    xlabel('U_{appl}/V');
    ylabel('Q_{ext}/C');
    % title('applied voltage over total charge');
    box on;set(gca,'Layer','top');
    saveas(gcf,lname,'bmp');
%     print(lname,'-dpdf','-noui','-bestfit');
%     savefig(sprintf('lissajous%s.fig',h));
    hold off; close(f);
    
    f = figure;
    hold on;
    x = time_volt(1:999)/1e-6;
    yyaxis left
    plot(time_volt(1:999)/1e-6,volt_gap(1:999),'k',time_volt(1:999)/1e-6,volt_appl(1:999),'k-.');
    xlabel('time in µs');
    ylabel('voltage in V');
    axis([min(x) max(x) -250 1250]);
    
    yyaxis right
    plot(time_volt(1:999)/1e-6,current_dis(1:999)*1000,'r');
    ylabel('current in mA');
    
    axis([min(x) max(x) -5 25]);
    zline = refline(0,0);
    zline.Color = 'k';
    zline.LineStyle = ':';
    
    legend('U_{gap}','U_{app}','I_{dis}');
    box on;set(gca,'Layer','top');
    % title('current/appl. & gap voltage via time');
    saveas(gcf,cname,'bmp');
%     print(cname,'-dpdf','-noui','-bestfit');
%     savefig(sprintf('currentdis%s.fig',h));
    hold off; close(f);
    
end

% Mache die Stark-Aufspaltung.

clear trsh Ind I inmax max nmax

[~, Ind] = max(abs(stark71_korr));
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
    plot(wavelength,starkshift-min(starkshift));
    x = wavelength;
    y = starkshift-min(starkshift);
    axis([min(x) max(x) min(y)-0.1*max(abs(y)) max(y)+0.1*max(abs(y))]);
    xlabel('wavelength in nm');
    ylabel('intensity, a.u.');
    box on;set(gca,'Layer','top');
    % title(sprintf('shift, 7.1in, %gµs, E=%g kV/cm', time, fieldstrength));
    saveas(gcf,'stark_shift71in','bmp');
%     print('stark_shift71in','-dpdf','-noui','-bestfit');
%     savefig(sprintf('stark_shift71in',h));
    hold off; close(f);

%% Daten.

save('2016-06-22.mat');

end