function stark_split2

warning off;
ppr_size = [14.6 11.4];

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

loc_dat = 'D:\documents\git\ag_praktikum_2016\2016-06-27\Daten\';
loc_main = 'D:\documents\git\ag_praktikum_2016\2016-06-27\';

% loc_dat= '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-22/Daten';
% loc_main = '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-22';

% Beschaffe mir die Dimensionen.

cd(loc_dat);
tmp = importdata('16Jun50001_RTO.dat');
tmp = tmp.data;
[n channel] = size(tmp);
time_volt = linspace(-0.8*1e-6,1.2*1e-6,n);
time_delta = time_volt(2)-time_volt(1);
cd(loc_main);

% Datei-Nummern und Wellenlängen.

m=2160;
data_nmb = 1:1:m;
loop_nmb = 1:1:15;
wvl_nmb = 1:1:36;
vert_nmb = 1:1:4;

raw_data = zeros(n,m);

charge = zeros(n,m);

wavelength = 491.8:0.02:492.5;
vertical = 6.8:0.1:7.1;

% Entnehmen der Datein.
% Es kommen zuerst die 10 loops zu einer Wellenlänge. Dann erneut 10 loops
% für die nächste. Macht 360 Datein, bis eine Höhe abgearbeitet ist.

cd(loc_dat);

for h=vert_nmb-1;
    
    for w=wvl_nmb-1;
        
        for l=loop_nmb;

            file_nmb = num2str(50000+l+15*w+540*h);
            nmb = l+15*w+540*h;
            disp(file_nmb);
            file = strcat('16Jun',file_nmb,'_RTO.dat');
            tmp = importdata(file);
            tmp = tmp.data;

            raw_data(:,nmb) = tmp(:,channel-1)+0.8;
            
            charge(:,nmb) = tmp(:,channel-2)+0.16;
        
        end
        
    end
    
end

cd(loc_main);
save('base_dat.mat');
% load('base_dat.mat');

% Kontrolle des Experimentes
    
    f = figure;
    hold on;
    meshc(time_volt/1e-6,linspace(1,2160,2160),charge');view(2);
    xlabel('time in µs');
    ylabel('vertical slit pos. in inch');
    c = colorbar;
    c.Label.String = 'voltage in V';
    title('charge via experiment duration');
    set(gca,'YTick',[270 810 1350 1890]);
    set(gca,'YTickLabel',{'6.8' '6.9' '7.0' '7.1'});
    
    set(c,'fontsize',12);
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'charge_full','pdf');
%     print('charge_full2','-dpdf','-noui','-bestfit');
    
%     savefig('charge_full');
    hold off; close(f);
    
% FFT abziehen der ersten 200 Punkte je Messung.

for i=1:2160  
    
    tmp = raw_data(1:500,i);
    tmpvec = [tmp' tmp' tmp' tmp' tmp'];
    tmpvec = tmpvec(1:n);
    tmpfft = fft(tmpvec);
    raw_fft(:,i) = fft(raw_data(:,i))-tmpfft';
    raw_fft(10:end-10,i) = 0;
    raw_data(:,i) = real(ifft(raw_fft(:,i)));
    
end

for i=1:36;
    
    stark_68in(:,i) = 1/15*sum(raw_data(:,(i-1)*15+1:(i-1)*15+15),2);
    stark_69in(:,i) = 1/15*sum(raw_data(:,(i-1)*15+540:(i-1)*15+555),2);
    stark_7in(:,i) = 1/15*sum(raw_data(:,(i-1)*15+1080:(i-1)*15+1095),2);
    stark_71in(:,i) = 1/15*sum(raw_data(:,(i-1)*15+1620:(i-1)*15+1635),2);
    
end

%% Bilder der Stark-Aufspaltungen.

% Wellenlängen-Offset. SGOLAY-Filt.
% Abziehen von Wellenlängen, bei denen keine Emission stattfindet.

    wvl68 = 1/8*sum(stark_68in(:,1:8),2);
    wvl69 = 1/8*sum(stark_69in(:,1:8),2);
    wvl7 = 1/8*sum(stark_7in(:,1:8),2);
    wvl71 = 1/8*sum(stark_71in(:,1:8),2);

    for i=1:36
        stark68_korr(:,i) = stark_68in(:,i)-wvl68;
        stark69_korr(:,i) = stark_69in(:,i)-wvl69;
        stark7_korr(:,i) = stark_7in(:,i)-wvl7;
        stark71_korr(:,i) = stark_71in(:,i)-wvl71;
    end
    
% SGOLAY-Filt.

    sgf_71in = sgolayfilt(stark71_korr,2,101);
    sgf_7in = sgolayfilt(stark7_korr,2,101);
    sgf_69in = sgolayfilt(stark69_korr,2,101);
    sgf_68in = sgolayfilt(stark68_korr,2,101); 
    
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
meshc(time_volt/1e-6,wavelength,sgf_69in');view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
title('vertical position 6.9 inch');
c = colorbar;
c.Label.String = 'intensity, a.u.';    

    set(c,'fontsize',12);
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_69in','pdf');
%     print('stark_69in2','-dpdf','-noui','-bestfit');

% savefig('stark_69in.fig');
hold off; close(f);

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,sgf_7in');view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
title('vertical position 7 inch');
c = colorbar;
c.Label.String = 'intensity, a.u.';   

    set(c,'fontsize',12);
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_7in','pdf');
%     print('stark_7in2','-dpdf','-noui','-bestfit');
    
% savefig('stark_7in.fig');
hold off; close(f);

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,sgf_71in');view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
title('vertical position 7.1 inch');c = colorbar;
c.Label.String = 'intensity, a.u.';
    
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'stark_71in','pdf');
%     print('stark_71in2','-dpdf','-noui','-bestfit');
    
% savefig('stark_71in.fig');
hold off; close(f);
    
% Schaue mir die Entladungscharakteristik an.

    file = strcat('16JunElectrics_RTO.dat');
    tmp = importdata(file);
    cd(loc_main);

    tmp = tmp.data;
    chrg_tmp(:) = tmp(:,3);
    volt_appl(:) = tmp(:,2);
    volt_diff(:) = 1/time_delta*diff(tmp(:,2));
    chrg_diff(:) = 1/time_delta*diff(tmp(:,3));
    rogowski(:) = tmp(:,5);
   
chrg_tmp = chrg_tmp*c_ext;
chrg_diff = c_ext*chrg_diff;
    
%Fit des C_tot und Berechnung des Entladungsstromes.

    xdata = volt_appl(1:150);
    ydata = chrg_tmp(1:150);
    
    P = polyfit(xdata,ydata,1);
    c_tot = P(1);
    
    c_par = c_tot-c_bd;
    
    volt_gap = volt_appl*(1-c_par/c_diel)-chrg_tmp./c_diel;
    current_dis = (1+c_gap/c_diel)*(chrg_diff-c_tot*volt_diff);

%Check.

    f = figure;
    hold on;
    plot(volt_appl,chrg_tmp);
    xlabel('U_{appl}/V');
    ylabel('Q_{ext}/C');
    title('applied voltage over totale charge');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'lissajous','bmp');
%     print('lissajous2','-dpdf','-noui','-bestfit');
%     savefig('lissajous.fig');
    hold off; close(f);
    
    f = figure;
    hold on;
    plot(time_volt(1:n-1)/1e-6,volt_gap(1:n-1),'r',time_volt(1:n-1)/1e-6,volt_appl(1:n-1),'r-.');
    xlabel('time in µs');
    yyaxis left
    ylabel('voltage in V');
    
    plot(time_volt(1:n-1)/1e-6,current_dis(1:n-1)*100000,'c');
    yyaxis right
    ylabel('current in mA');
    
    legend('U_{gap}','U_{app}','I_{dis}');
    title('current/appl. & gap voltage via time');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'currentdis','bmp');
%     print('currentdis2','-dpdf','-noui','-bestfit');
%     savefig('currentdis.fig');
    hold off; close(f);

%% Daten.

save('2016-06-27.mat');

end