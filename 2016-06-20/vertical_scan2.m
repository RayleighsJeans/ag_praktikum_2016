function vertical_scan2;

warning off;

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

% loc_dat= '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-15/Daten';
% load '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-15/base_data.mat';

loc_dat = 'D:\documents\git\ag_praktikum_2016\2016-06-20\Daten';
loc_main = 'D:\documents\git\ag_praktikum_2016\2016-06-20\';
load 'D:\documents\git\ag_praktikum_2016\2016-06-20\base_dat.mat';

% Datei-Nummern und Wellenlängen. Dimensionen und Initialisierung der späteren Messwerte-Matrizen.

m = 1450;
a=30000;

    dat_587 = 1+a:50:m+a;
    dat_667 = 11+a:50:m+a;
    dat_690 = 21+a:50:m+a;
    dat_706 = 31+a:50:m+a;
    dat_728 = 41+a:50:m+a;
    
data_nmb = [dat_587' dat_667' dat_690' dat_706' dat_728'];
    
    current = zeros(1999,m);
    chrg = zeros(2000,m);

wavelength = [587.65 667.96 690.0 706.66 728.31];

cd(loc_dat);
tmp = importdata('16JunElectrics_RTO.dat');
tmp = tmp.data;
time_volt = tmp(:,1);
time_delta = time_volt(2)-time_volt(1);
cd(loc_main);

% Vertikale Position.

vertical_pos = 5.8:0.05:7.2;
vertical_mm = -0.2:0.117:3.1;

% Entnehme die Daten aus den Files.

cd(loc_dat);

for i=1:29
    
    for k=1:10
        
    itr = data_nmb(i,:);
        
%     Entnehme Daten aus den Datein. Hier für 690nm.
    nmb = num2str(itr(3)-1+k);
    file = strcat('16Jun',nmb,'_RTO.dat');
    disp(file)
    C = importdata(file);
    tmp_nm690 = C.data;

        nm690(:,(i-1)*10+k) = tmp_nm690(:,4);
        current(:,itr(3)-1+k) = 1/time_delta*diff(tmp_nm690(:,3));
        chrg(:,itr(3)-1+k) = tmp_nm690(:,3);
        
%     Hier für 587.65nm.
    nmb = num2str(itr(1)-1+k);
    file = strcat('16Jun',nmb,'_RTO.dat');
    disp(file)
    A = importdata(file);
    tmp_nm587 = A.data;
    
        nm587(:,(i-1)*10+k) = tmp_nm587(:,4);
        current(:,itr(1)-1+k) = 1/time_delta*diff(tmp_nm587(:,3));
        chrg(:,itr(1)-1+k) = tmp_nm587(:,3);
    
%     Hier für 667.96nm.´
    nmb = num2str(itr(2)-1+k);
    file = strcat('16Jun',nmb,'_RTO.dat');
    disp(file)
    B = importdata(file);
    tmp_nm667 = B.data;

        nm667(:,(i-1)*10+k) = tmp_nm667(:,4);
        current(:,itr(2)-1+k) = 1/time_delta*diff(tmp_nm667(:,3));
        chrg(:,itr(2)-1+k) = tmp_nm667(:,3);

%     Hier für 706.66nm.
    nmb = num2str(itr(4)-1+k);
    file = strcat('16Jun',nmb,'_RTO.dat');
    disp(file)
    D = importdata(file);
    tmp_nm706 = D.data;

    nm706(:,(i-1)*10+k) = tmp_nm706(:,4);
    current(:,itr(4)-1+k) = 1/time_delta*diff(tmp_nm706(:,3));
        chrg(:,itr(4)-1+k) = tmp_nm706(:,3);

%     Hier für 728.31nm.

    nmb = num2str(itr(5)-1+k);
    file = strcat('16Jun',nmb,'_RTO.dat');
    disp(file)
    E = importdata(file);
    tmp_nm728 = E.data;

        nm728(:,(i-1)*10+k) = tmp_nm728(:,4);
        current(:,itr(5)-1+k) = 1/time_delta*diff(tmp_nm728(:,3));
        chrg(:,itr(5)-1+k) = tmp_nm728(:,3);

    end
    
end

% Zurück.

% cd(loc_main);

        %Mache Off-Set aus der 690nm-Messung. Nur ein mal.
        for i=1:29
            offset(:,i) = 1/10*sum(nm690(:,(i-1)*10+1:(i-1)*10+10),2);
        end
        
% Offset abziehen.

for i=1:29
    for k=1:10 
     nm587(:,(i-1)*10+k) = nm587(:,(i-1)*10+k)-offset(:,i)-max(max(nm587(:,(i-1)*10+k)-offset(:,i)));
     nm667(:,(i-1)*10+k) = nm667(:,(i-1)*10+k)-offset(:,i)-max(max(nm667(:,(i-1)*10+k)-offset(:,i)));
     nm706(:,(i-1)*10+k) = nm706(:,(i-1)*10+k)-offset(:,i)-max(max(nm706(:,(i-1)*10+k)-offset(:,i)));
     nm728(:,(i-1)*10+k) = nm728(:,(i-1)*10+k)-offset(:,i)-max(max(nm728(:,(i-1)*10+k)-offset(:,i)));
    end
end

    %Mittelungen.
    
    for i=1:29

        mean587(:,i) = 1/10*sum(nm587(:,(i-1)*10+1:(i-1)*10+10),2);
        mean667(:,i) = 1/10*sum(nm667(:,(i-1)*10+1:(i-1)*10+10),2);
        mean706(:,i) = 1/10*sum(nm706(:,(i-1)*10+1:(i-1)*10+10),2);
        mean728(:,i) = 1/10*sum(nm728(:,(i-1)*10+1:(i-1)*10+10),2);

    end

cd(loc_main);
save base_dat.mat
% load base_dat.mat

% Bilder

f = figure;hold on;
meshc(time_volt/1e-6,vertical_pos,real(log10(mean587))');view(2);
% title('vertical emission profile at 587.65 nm');
c = colorbar;
c.Label.String = 'log_{10} of intensity, a.u.';
ylabel('vertical pos. in inch');
ax = gca; ax.YTickLabel = {'anode','6','6.2','6.4','6.6','6.8','7','kathode'};
ax.YTickLabelRotation = 90;
xlabel('time in µs');
box on;set(gca,'Layer','top');
% savefig('587nm.fig');
%     saveas(gcf,'587nm','pdf');
    print('587nm','-dpdf','-noui','-bestfit');

hold off;close(f);

f = figure;hold on;
meshc(time_volt/1e-6,vertical_pos,real(log10(mean667))');view(2);
% title('vertical emission profile for 667.96 nm');
c = colorbar;
c.Label.String = 'log_{10} of intensity, a.u.';
ylabel('vertical pos. in inch');
ax = gca; ax.YTickLabel = {'anode','6','6.2','6.4','6.6','6.8','7','kathode'};
ax.YTickLabelRotation = 90;
xlabel('time in µs');
box on;set(gca,'Layer','top');   
% savefig('667nm.fig');
%     saveas(gcf,'667nm','pdf');
    print('667nm','-dpdf','-noui','-bestfit');

hold off;close(f);

f = figure;hold on;
meshc(time_volt/1e-6,vertical_pos,real(log10(mean706))');view(2);
% title('vertical emission profile for 706.66 nm');
c = colorbar;
c.Label.String = 'log_{10} of intensity, a.u.';
ylabel('vertical pos. in inch');
ax = gca; ax.YTickLabel = {'anode','6','6.2','6.4','6.6','6.8','7','kathode'};
ax.YTickLabelRotation = 90;
xlabel('time in µs');
box on;set(gca,'Layer','top'); 
% savefig('706nm.fig');

%     saveas(gcf,'706nm','pdf');
    print('706nm','-dpdf','-noui','-bestfit');

hold off;close(f);

f = figure;hold on;
meshc(time_volt/1e-6,vertical_pos,real(log10(mean728))');view(2);
% title('vertical emission profile for 728.31 nm');
c = colorbar;
c.Label.String = 'log_{10} of intensity, a.u.';
ylabel('vertical pos. in inch');
ax = gca; ax.YTickLabel = {'anode','6','6.2','6.4','6.6','6.8','7','kathode'};
ax.YTickLabelRotation = 90;
xlabel('time in µs');
box on;set(gca,'Layer','top'); 
% savefig('728nm.fig');

%     saveas(gcf,'728nm','pdf');
    print('728nm','-dpdf','-noui','-bestfit');

hold off;close(f);

% Schaue mir die Entladungscharakteristik an.

cd(loc_dat);
tmp = importdata('16JunElectrics_RTO.dat');
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

    f = figure;hold on; box on;
    plot(volt_appl,chrg_tmp);
    x = volt_appl;
    y = chrg_tmp;
    axis([min(x) max(x) min(y)-0.1*max(abs(y)) max(y)+0.1*max(abs(y))]);
    xlabel('U_{appl}/V');
    ylabel('Q_{ext}/C');
    box on;set(gca,'Layer','top'); 
    % title('applied voltage over total charge');
%     savefig('lissajous.fig');
    
        saveas(gcf,'lissajous','bmp');
%         print('lissajous2','-dpdf','-noui','-bestfit');
    
    hold off;close(f);

    
    f = figure;
    hold on;
    xlabel('time in µs');
    yyaxis left
    plot(time_volt(1:1999)/1e-6,volt_gap(1:1999),'k-',time_volt(1:1999)/1e-6,volt_appl(1:1999),'k-.');
    x = time_volt/1e-6;
    ylabel('voltage in V');
    axis([min(x) max(x) -250 1250]);
    
    yyaxis right
    plot(time_volt(1:1999)/1e-6,current_dis(1:1999)*1000,'r-');
    ylabel('current in mA');
    x = time_volt/1e-6;
    
    legend('U_{gap}','U_{app}','I_{dis}');
    box on;set(gca,'Layer','top'); 
    zline = refline(0,0);
    zline.Color = 'k';
    zline.LineStyle = ':';
    % title('current/appl. & gap voltage via time');
%     savefig('current_dis.fig');
    
        saveas(gcf,'current_dis','bmp');
%         print('current_dis2','-dpdf','-noui','-bestfit');
    
    hold off;close(f);
    
% Linienverhältnisse.

    tmp706 = mean706;
    tmp587 = mean587;
    tmp667 = mean667;
    tmp728 = mean728;

    schwell706 = 0.05*max(max(abs(tmp706)));
    schwell587 = 0.05*max(max(abs(tmp587)));
    schwell667 = 0.05*max(max(abs(tmp667)));
    schwell728 = 0.05*max(max(abs(tmp728)));
    
    for i=1:2000
        for j=1:29
       
            if (abs(tmp706(i,j))<=schwell706)
                tmp706(i,j)=0;
            end
            
%             if (tmp587(i,j)<=schwell587)
%                 tmp587(i,j)=0;
%             end
            
            if (abs(tmp667(i,j))<=schwell667)
                tmp667(i,j)=0;
            end
            
%             if (tmp728(i,j)<=schwell728)
%                 tmp728(i,j)=0;
%             end
                        
        end
    end
    
    % Numerischer Umweg.
    
        tmp667 = tmp667*1e4;
        tmp728 = tmp728*1e4;
        
        tmp706 = tmp706*1e3;
        tmp587 = tmp587*1e4;
    
    tmp_ratio706 = tmp706./tmp587;
    tmp_ratio667 = tmp667./tmp728;
    
    sgf_ratio706 = sgolayfilt(tmp_ratio706,2,51);
    sgf_ratio667 = sgolayfilt(tmp_ratio667,2,51);
    
    field667 = tmp_ratio667.^2*45.07-20.18*tmp_ratio667-19.98*tmp_ratio667.^3+3.369*tmp_ratio667.^4+2.224;
    sgf_field667 = sgolayfilt(field667,2,51);
    field706 = tmp_ratio706.^2*45.07-20.18*tmp_ratio706-19.98*tmp_ratio706.^3+3.369*tmp_ratio706.^4+2.224;
    sgf_field706 = sgolayfilt(field706,2,51);
    
    f = figure; hold on;
    meshc(time_volt/1e-6,vertical_pos,sgf_ratio706');view(2);
    ylabel('vertical slit pos. in inch');
    ax = gca; ax.YTickLabel = {'anode','6','6.2','6.4','6.6','6.8','7','kathode'};
    ax.YTickLabelRotation = 90;
    c = colorbar;
    c.Label.String = 'line ratio, a.u.';
    box on;set(gca,'Layer','top');
    xlabel('time in µs');
    % title('line ratio of He lines at 706 nm and 587 nm');
%     savefig('lineratio706.fig');
    
%         saveas(gcf,'lineratio706','pdf');
        print('lineratio706','-dpdf','-noui','-bestfit');
    
    hold off;close(f);
    
    f = figure; hold on;
    meshc(time_volt/1e-6,vertical_pos,sgf_ratio667');view(2);
    ylabel('vertical slit pos. in inch');
    ax = gca; ax.YTickLabel = {'anode','6','6.2','6.4','6.6','6.8','7','kathode'};
    ax.YTickLabelRotation = 90;
    c = colorbar;
    c.Label.String = 'line ratio, a.u.';
    xlabel('time in µs');
    % title('line ratio of He lines at 667 nm and 728 nm');
    box on;set(gca,'Layer','top'); 
%     savefig('lineratio667.fig');
    
%         saveas(gcf,'lineratio667','pdf');
        print('lineratio667','-dpdf','-noui','-bestfit');
    
    hold off;close(f);

    f = figure; hold on;
    meshc(time_volt/1e-6,vertical_pos,sgf_field667');view(2);
    ylabel('vertical slit pos. in inch');
    ax = gca; ax.YTickLabel = {'anode','6','6.2','6.4','6.6','6.8','7','kathode'};
    ax.YTickLabelRotation = 90;
    c = colorbar;
    c.Label.String = 'el. field in kV/cm';
    xlabel('time in µs');
    % title('el. field strength from ratio 667nm/728nm');
    box on;set(gca,'Layer','top'); 
%     savefig('lineratio667.fig');
    
%         saveas(gcf,'elfield667','pdf');
        print('elfield667','-dpdf','-noui','-bestfit');
    
    hold off;close(f);

    f = figure; hold on;
    meshc(time_volt/1e-6,vertical_pos,sgf_field706');view(2);
    ylabel('vertical slit pos. in inch');
    ax = gca; ax.YTickLabel = {'anode','6','6.2','6.4','6.6','6.8','7','kathode'};
    ax.YTickLabelRotation = 90;
    c = colorbar;
    c.Label.String = 'el. field in kV/cm';
    xlabel('time in µs');
    % title('el. field strength from ratio 706nm/587nm');
    box on;set(gca,'Layer','top'); 
%     savefig('lineratio706.fig');
    
%         saveas(gcf,'elfield706','pdf');
        print('elfield706','-dpdf','-noui','-bestfit');
    
    hold off;close(f);
    
% Daten.
save('2016-06-20.mat');

end