function stark_split

set(0, 'DefaultTextInterpreter', 'latex')
warning off;
ppr_size = [14.6 11.4];

%% Initialisierung

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


%% Dateipfade der Messdaten.

loc_dat = 'D:\documents\git\ag_praktikum_2016\2016-06-22\Daten';
loc_main = 'D:\documents\git\ag_praktikum_2016\2016-06-22\';

% loc_dat= '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-22/Daten';
% loc_main = '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-22';

%% Beschaffe mir die Dimensionen.

cd(loc_dat);
tmp = importdata('16Jun35001_RTO.dat');
tmp = tmp.data;
[n channel] = size(tmp);
time_volt = tmp(:,1);
time_delta = time_volt(2)-time_volt(1);
cd(loc_main);

%% Datei-Nummern und Wellenlängen.

m=1440;
data_nmb = 1:1:m;
loop_nmb = 1:1:10;
wvl_nmb = 1:1:36;
vert_nmb = 1:1:4;

raw_data = zeros(n,m);

rogowski = zeros(n,m);

wavelength = 491.8:0.02:492.5;
vertical = 6.8:0.1:7.1;

%% Entnehmen der Datein.
%Es kommen zuerst die 10 loops zu einer Wellenlänge. Dann erneut 10 loops
%für die nächste. Macht 360 Datein, bis eine Höhe abgearbeitet ist.

cd(loc_dat);

for h=vert_nmb-1;
    
    for w=wvl_nmb-1;
        
        for l=loop_nmb;

            file_nmb = num2str(35000+l+10*w+360*h);
            nmb = l+10*w+360*h;
            disp(file_nmb);
            file = strcat('16Jun',file_nmb,'_RTO.dat');
            tmp = importdata(file);
            tmp = tmp.data;

            raw_data(:,nmb) = tmp(:,4)+0.16;
            
            rogowski(:,nmb) = tmp(:,5);
        
        end
        
    end
    
end

cd(loc_main);
save('base_dat.mat');

% Skippe die Schleifen und Initialiserung, lade Basis-Datensatz.
% load base_dat.mat

% Kontrolle des Experimentes
    
    f = figure;
    hold on;
    meshc(linspace(1,1440,1440),time_volt,rogowski+0.160);view(2);
    ylabel('time in µs');
    xlabel('vertical slit pos. in inch');
    c = colorbar;
    c.Label.String = 'voltage in V';
    title('rogowski coil via experiment duration');
    set(gca,'XTick',[180 540 900 1260]);
    set(gca,'XTickLabel',{'6.8' '6.9' '7.0' '7.1'});
    
    set(c,'fontsize',12);
    set(gca,'FontSize', 12,'Fontname','L M Roman12');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'rogowski_full','pdf');
    print('rogowski_full2','-dpdf','-noui','-bestfit');
    
    savefig('rogowski_full');
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
meshc(wavelength,time_volt/1e-6,sgf_68in);view(2);
xlabel('wavelength in nm');
ylabel('time in µs');
title('vertical position 6.8 inch');
c = colorbar;
c.Label.String = 'intensity, a.u.'; 

    set(c,'fontsize',12);
    set(gca,'FontSize', 12,'Fontname','L M Roman12');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_68in','pdf');
    print('stark_68in2','-dpdf','-noui','-bestfit');
 
savefig('stark_68in.fig');
hold off; close(f);

f = figure;hold on;
meshc(wavelength,time_volt/1e-6,sgf_69in);view(2);
xlabel('wavelength in nm');
ylabel('time in µs');
title('vertical position 6.9 inch');
c = colorbar;
c.Label.String = 'intensity, a.u.';    

    set(c,'fontsize',12);
    set(gca,'FontSize', 12,'Fontname','L M Roman12');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_69in','pdf');
    print('stark_69in2','-dpdf','-noui','-bestfit');
 
savefig('stark_69in.fig');
hold off; close(f);

f = figure;hold on;
meshc(wavelength,time_volt/1e-6,sgf_7in);view(2);
xlabel('wavelength in nm');
ylabel('time in µs');
title('vertical position 7 inch'); 
c = colorbar;
c.Label.String = 'intensity, a.u.';   

    set(c,'fontsize',12);
    set(gca,'FontSize', 12,'Fontname','L M Roman12');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_7in','pdf');
    print('stark_7in2','-dpdf','-noui','-bestfit');
 
savefig('stark_7in.fig');
hold off; close(f);

f = figure;hold on;
meshc(wavelength,time_volt/1e-6,sgf_71in);view(2);
xlabel('wavelength in nm');
ylabel('time in µs');
title('vertical position 7.1 inch');
c = colorbar;
c.Label.String = 'intensity, a.u.';
    
    set(gca,'FontSize', 12,'Fontname','L M Roman12');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_71in','pdf');
    print('stark_71in2','-dpdf','-noui','-bestfit');
 
savefig('stark_71in.fig');
hold off; close(f);

% Schaue mir die Entladungscharakteristik an.
% load 2016-06-22.mat

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

    lname1 = strcat('lissajous',h,'.pdf');
    lname2 = strcat('lissajous2',h,'.pdf');
    
    cname1 = strcat('currentdis',h,'.pdf');
    cname2 = strcat('currentdis2',h,'.pdf');

    f = figure;hold on;
    plot(volt_appl,chrg_tmp);
    xlabel('U_{appl}/V');
    ylabel('Q_{ext}/C');
    title('applied voltage over total charge');
    set(gca,'FontSize', 12,'Fontname','L M Roman12');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,lname1,'pdf');
    print(lname2,'-dpdf','-noui','-bestfit');
    savefig(sprintf('lissajous%s.fig',h));
    hold off; close(f);
    
    f = figure;
    hold on;
    plot(time_volt(1:999)/1e-6,volt_gap(1:999),'r',time_volt(1:999)/1e-6,volt_appl(1:999),'r-.');
    xlabel('time in µs');
    yyaxis left
    ylabel('voltage in V');
    
    plot(time_volt(1:999)/1e-6,current_dis(1:999)*100000,'c');
    yyaxis right
    ylabel('current in 10^{-5} A');
    
    legend('U_{gap}','U_{app}','I_{dis}');
    title('current/appl. & gap voltage via time');
    set(gca,'FontSize', 12,'Fontname','L M Roman12');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,cname1,'pdf');
    print(cname2,'-dpdf','-noui','-bestfit');
    savefig(sprintf('currentdis%s.fig',h));
    hold off; close(f);
    
end

%% Daten.

save('2016-06-22.mat');

end