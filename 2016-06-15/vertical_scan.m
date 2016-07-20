function vertical_scan

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

%loc_dat= '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-15/Daten';
%load '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-15/offset_data.mat';

loc_dat = 'D:\documents\git\ag_praktikum_2016\2016-06-15\Daten\';
load 'D:\documents\git\ag_praktikum_2016\2016-06-15\offset_data.mat';
loc_main = 'D:\documents\git\ag_praktikum_2016\2016-06-15\';


%% Datei-Nummern und Wellenlängen.

    dat_587 = 1:25:725;
    dat_667 = 6:25:725;
    dat_690 = 11:25:725;
    dat_706 = 16:25:725;
    dat_728 = 21:25:725;
    
data_nmb = [dat_587' dat_667' dat_690' dat_706' dat_728'];
    
chrg_current = zeros(1999,725);
chrg = zeros(2000,725);
volt_diff = zeros(1999,725);
volt_appl = zeros(2000,725);
control = zeros(2000,725);

wavelength = [587.65 667.96 690.0 706.66 728.31];

cd(loc_dat);
tmp = importdata('16Jun16001_RTO.dat');
tmp = tmp.data;
time_volt = tmp(:,1);
time_delta = time_volt(2)-time_volt(1);

%% Vertikaler Verfahr-Weg.

vertical_pos = 5.8:0.05:7.2;
vertical_mm = -0.2:0.117:3.1;

cd(loc_dat);

%% Entnehmen der Datein.

for k=1:5
    
    for i=1:29
        
    itr = data_nmb(i,:);
    
    nmb = num2str(itr(3)-1+k,'%03d');
    file = strcat('16Jun16',nmb,'_RTO.dat');
    disp(file)
    C = importdata(file);
    tmp_nm690 = C.data;

    nm690(:,(i-1)*5+k) = tmp_nm690(:,4);
    
        volt_appl(:,itr(3)-1+k) = tmp_nm690(:,2);
        volt_diff(:,itr(3)-1+k) = 1/time_delta*diff(tmp_nm690(:,2));
    
        control(:,itr(3)-1+k) = tmp_nm690(:,4);
        chrg_current(:,itr(3)-1+k) = 1/time_delta*diff(tmp_nm690(:,3));
        chrg(:,itr(3)-1+k) = tmp_nm690(:,3);
        
    %587nm
    nmb = num2str(itr(1)-1+k,'%03d');
    file = strcat('16Jun16',nmb,'_RTO.dat');
    disp(file)
    A = importdata(file);
    tmp_nm587 = A.data;
    
    nm587(:,(i-1)*5+k) = tmp_nm587(:,4)-offset(:,i)-max(max(tmp_nm587(:,4)-offset(:,i)));
    
        volt_appl(:,itr(1)-1+k) = tmp_nm587(:,2);
        volt_diff(:,itr(1)-1+k) = 1/time_delta*diff(tmp_nm587(:,2));       
    
        control(:,itr(1)-1+k) = tmp_nm587(:,4);
        chrg_current(:,itr(1)-1+k) = 1/time_delta*diff(tmp_nm587(:,3));
        chrg(:,itr(1)-1+k) = tmp_nm587(:,3);
    
    %667nm
    nmb = num2str(itr(2)-1+k,'%03d');
    file = strcat('16Jun16',nmb,'_RTO.dat');
    disp(file)
    B = importdata(file);
    tmp_nm667 = B.data;

    nm667(:,(i-1)*5+k) = tmp_nm667(:,4)-offset(:,i)-max(max(tmp_nm667(:,4)-offset(:,i)));
    
        volt_appl(:,itr(2)-1+k) = tmp_nm667(:,2);
        volt_diff(:,itr(2)-1+k) = 1/time_delta*diff(tmp_nm667(:,2));
        
        control(:,itr(2)-1+k) = tmp_nm667(:,4);
        chrg_current(:,itr(2)-1+k) = 1/time_delta*diff(tmp_nm667(:,3));
        chrg(:,itr(2)-1+k) = tmp_nm667(:,3);

    %706nm
    nmb = num2str(itr(4)-1+k,'%03d');
    file = strcat('16Jun16',nmb,'_RTO.dat');
    disp(file)
    D = importdata(file);
    tmp_nm706 = D.data;

    nm706(:,(i-1)*5+k) = tmp_nm706(:,4)-offset(:,i)-max(max(tmp_nm706(:,4)-offset(:,i)));
    
        volt_appl(:,itr(4)-1+k) = tmp_nm706(:,2);
        volt_diff(:,itr(4)-1+k) = 1/time_delta*diff(tmp_nm706(:,2));
    
        control(:,itr(4)-1+k) = tmp_nm706(:,4);
        chrg_current(:,itr(4)-1+k) = 1/time_delta*diff(tmp_nm706(:,3));
        chrg(:,itr(4)-1+k) = tmp_nm706(:,3);

    %728nm
    nmb = num2str(itr(5)-1+k,'%03d');
    file = strcat('16Jun16',nmb,'_RTO.dat');
    disp(file)
    E = importdata(file);
    tmp_nm728 = E.data;

    nm728(:,(i-1)*5+k) = tmp_nm728(:,4)-offset(:,i)-max(max(tmp_nm728(:,4)-offset(:,i)));
    
        volt_appl(:,itr(5)-1+k) = tmp_nm728(:,2);
        volt_diff(:,itr(5)-1+k) = 1/time_delta*diff(tmp_nm728(:,2));
    
        control(:,itr(5)-1+k) = tmp_nm728(:,4);
        chrg_current(:,itr(5)-1+k) = 1/time_delta*diff(tmp_nm728(:,3));
        chrg(:,itr(5)-1+k) = tmp_nm728(:,3);
    
    end
    
end

%% Berechnen der sauberen Spektren.

for i=1:29
    
    mean587(:,i) = 1/5*sum(nm587(:,(i-1)*5+1:(i-1)*5+5),2);
    mean667(:,i) = 1/5*sum(nm667(:,(i-1)*5+1:(i-1)*5+5),2);
    mean706(:,i) = 1/5*sum(nm706(:,(i-1)*5+1:(i-1)*5+5),2);
    mean728(:,i) = 1/5*sum(nm728(:,(i-1)*5+1:(i-1)*5+5),2);
    
end

cd(loc_main);

%% Bilder der vertikal aufgelösten Wellenlängen.

f = figure;
hold on;
meshc(time_volt/1e-6,vertical_pos,real(log10(mean587))');view(2);
title('vertical emission profile at 587.65 nm');
c = colorbar;
c.Label.String = 'log_{10} of intensity, a.u.';
ylabel('vertical slit pos. in inch');
xlabel('time in µs');
set(c,'fontsize',12);
% savefig('korr587nm.fig');

    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'korr587nm','pdf');
%     print('korr587nm2','-dpdf','-noui','-bestfit');

hold off;close(f);


f = figure;
hold on;
meshc(time_volt/1e-6,vertical_pos,real(log10(mean667))');view(2);
title('vertical emission profile at 667.96 nm');
c = colorbar;
c.Label.String = 'log_{10} of intensity, a.u.';
ylabel('vertical slit pos. in inch');
xlabel('time in µs');
set(c,'fontsize',12);
% savefig('korr667nm.fig');

    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'korr667nm','pdf');
%     print('korr667nm2','-dpdf','-noui','-bestfit');

hold off;close(f);

f = figure;
hold on;
meshc(time_volt/1e-6,vertical_pos,real(log10(offset))');view(2);
title('vertical emission profile at 690 nm');
c = colorbar;
c.Label.String = 'log_{10} of intensity, a.u.';
ylabel('vertical slit pos. in inch');
xlabel('time in µs');
set(c,'fontsize',12);
% savefig('korr690nm.fig');

    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'korr690nm','pdf');
%     print('korr690nm2','-dpdf','-noui','-bestfit');

hold off;close(f);

f = figure;
hold on;
meshc(time_volt/1e-6,vertical_pos,real(log10(mean706))');view(2);
title('vertical emission profile at 706.66 nm');
c = colorbar;
c.Label.String = 'log_{10} of intensity, a.u.';
ylabel('vertical slit pos. in inch');
xlabel('time in µs');
set(c,'fontsize',12);
    
% savefig('korr706nm.fig');

    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'korr706nm','pdf');
%     print('korr706nm2','-dpdf','-noui','-bestfit');

hold off;close(f);

f = figure;
hold on;
meshc(time_volt/1e-6,vertical_pos,real(log10(mean728))');view(2);
title('vertical emission profile at 728.34 nm');
c = colorbar;
c.Label.String = 'log_{10} of intensity, a.u.';
ylabel('vertical slit pos. in inch');
xlabel('time in µs');
set(c,'fontsize',12);
    
% savefig('korr728nm.fig');

    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'korr728nm','pdf');
%     print('korr728nm2','-dpdf','-noui','-bestfit');

hold off;close(f);

chrg = chrg*c_ext;
chrg_current = c_ext*chrg_current;

%% Fit des C_tot und Berechnung des Entladungsstromes.

    xdata = volt_appl(1:150,1);
    ydata = chrg(1:150,1);
    
    P = polyfit(xdata,ydata,1);
    c_tot = P(1);
    
    c_par = c_tot-c_bd;
    
    volt_gap = volt_appl*(1-c_par/c_diel)-chrg./c_diel;
    current_dis = (1+c_gap/c_diel)*(chrg_current-c_tot*volt_diff);

%load 2016-06-15

    f = figure;
    hold on;
    plot(volt_appl(:,1),chrg(:,1));
    xlabel('U_{appl}/V');
    ylabel('Q_{ext}/C');
    title('applied voltage over totale charge');
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'lissajous','bmp');
%     print('lissajous2','-dpdf','-noui','-bestfit');
%     savefig('lissajous.fig');
    hold off;close(f);
    

    f = figure;
    hold on;
    plot(time_volt(1:1999)/1e-6,volt_gap(1:1999,1),'r',time_volt(1:1999)/1e-6,volt_appl(1:1999,1),'r-.');
    xlabel('time in µs');
    yyaxis left
    ylabel('voltage in V');
    
    plot(time_volt(1:1999)/1e-6,current_dis(1:1999,1)*100000,'c');
    yyaxis right
    ylabel('current in mA');
    
    legend('U_{gap}','U_{app}','I_{dis}');
    title('current/appl. & gap voltage via time');
    
    set(gcf,'PaperSize',ppr_size);
    saveas(gcf,'current_dis','bmp');
%     print('current_dis2','-dpdf','-noui','-bestfit');
    
%     savefig('current_dis.fig');
    
    hold off;close(f);
    
%% Linienverhältnisse.
load 2016-06-15.mat

    tmp1 = mean706./mean587;
    tmp2 = mean667./mean728;

    schwell1 = 0.3*max(max(tmp1-min(min(tmp1))));
    schwell2 = 0.25*max(max(tmp2-min(min(tmp2))));
    
    tmp1 = tmp1-min(min(tmp1));
    tmp2 = tmp2-min(min(tmp2));
    
    for i=1:2000
        for j=1:29
       
            if (tmp1(i,j)<=schwell1)
                tmp1(i,j)=0;
            end
            
            if (tmp2(i,j)<=schwell2)
                tmp2(i,j)=0;
            end
            
        end
    end

    tmp1 = sgolayfilt(tmp1,2,51);
    tmp2 = sgolayfilt(tmp2,2,51);
    
    f = figure; hold on;
    meshc(time_volt/1e-6,vertical_pos,tmp1');view(2);
    ylabel('vertical slit pos. in inch');
    c = colorbar;
    c.Label.String = 'line ratio, a.u.';
    xlabel('time in µs');
    
    set(c,'fontsize',12);
    title('line ratio of He lines at 706 nm and 587 nm');
    
        set(gcf,'PaperSize',ppr_size);
        saveas(gcf,'lineratio706','pdf');
%         print('lineratio7062','-dpdf','-noui','-bestfit');
        
%     savefig('lineratio706.fig');
    hold off;close(f);
    
    f = figure; hold on;
    meshc(time_volt/1e-6,vertical_pos,tmp2');view(2);
    ylabel('vertical slit pos. in inch');
    c = colorbar;
    c.Label.String = 'line ratio, a.u.';
    xlabel('time in µs');
    
    set(c,'fontsize',12);
    title('line ratio of He lines at 667 nm and 728 nm');
    
        set(gcf,'PaperSize',ppr_size);
        saveas(gcf,'lineratio728','pdf');
%         print('lineratio7282','-dpdf','-noui','-bestfit');
    
%     savefig('lineratio728.fig');
    hold off;close(f);
    
%% Daten.

save('2016-06-15.mat');

end