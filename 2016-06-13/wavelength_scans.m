function wavelength_scans;

warning off;

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

%Dateipfade der Messdaten

loc_main = 'D:\documents\git\ag_praktikum_2016\2016-06-13\';
loc_a = 'D:\documents\git\ag_praktikum_2016\2016-06-13\587nm';
loc_b = 'D:\documents\git\ag_praktikum_2016\2016-06-13\667nm';
loc_c = 'D:\documents\git\ag_praktikum_2016\2016-06-13\706nm';
loc_d = 'D:\documents\git\ag_praktikum_2016\2016-06-13\728nm';
loc_e = 'D:\documents\git\ag_praktikum_2016\2016-06-13\spectrum';

%Messbereiche abgesteckt

scan1 = linspace(587,588,51);
scan2 = linspace(667.3,668.3,51);
scan3 = linspace(706,707,51);
scan4 = linspace(727.7,728.7,51);
scan_spec = linspace(585,730,291);

cd(loc_a);
tmp = importdata('16Jun08001_RTO.dat');
tmp = tmp.data;
time_volt = tmp(:,1);
time_delta = time_volt(2)-time_volt(1);
cd(loc_main);

    % Initialisierung der Messwerte-Matrizen
    
    nm587 = zeros(1000,51);
    nm667 = zeros(1000,51);
    nm706= zeros(1000,51);
    nm728 = zeros(1000,51);
    spectrum = zeros(1000,291);
    volt_appl = zeros(1000,495);
    volt_diff = zeros(1000,495);
    chrg = zeros(1000,495);
    chrg_current = zeros(1000,495);
    
%% Hole mir die Messdaten aus den RTO-files. Nehme den Photomultiplier-Strom und mache direkt die Ableitung vom Entladungstrom.
    
for i=1:51;
    
    cd(loc_a);
    nmb = num2str(i,'%02d');
    file_a = strcat('16Jun080',nmb,'_RTO.dat');
    disp(file_a);
    A = importdata(file_a);
    tmp_nm587 = A.data;

    %Beispielhaft Daten entnommen. Schreibe in Ergebnis-Matrix/-Vektor.
    
    nm587(:,i) = tmp_nm587(:,4);
    
        chrg_current(1000,i) = tmp_nm587(1000,3);
        chrg(:,i) = tmp_nm587(:,3);
        chrg_current(1:999,i) = 1/time_delta*diff(tmp_nm587(:,3));

        volt_appl(:,i) = tmp_nm587(:,2);
        volt_diff(1000,i) = tmp_nm587(1000,2);
        volt_diff(1:999,i) = 1/time_delta*diff(tmp_nm587(:,2));
    
    cd(loc_b);
    file_b = strcat('16Jun081',nmb,'_RTO.dat');
    disp(file_b);
    B = importdata(file_b);
    tmp_nm667 = B.data;
    
    nm667(:,i) = tmp_nm667(:,4);
    
        chrg_current(1000,i+51) = tmp_nm667(1000,3);
        chrg(:,i+51) = tmp_nm667(:,3);
        chrg_current(1:999,i+51) = 1/time_delta*diff(tmp_nm667(:,3));

        volt_appl(:,i+51) = tmp_nm667(:,2);
        volt_diff(1000,i+51) = tmp_nm667(1000,2);
        volt_diff(1:999,i+51) = 1/time_delta*diff(tmp_nm667(:,2));
    
    cd(loc_c);
    file_c = strcat('16Jun082',nmb,'_RTO.dat');
    disp(file_c);
    C = importdata(file_c);
    tmp_nm706 = C.data;
    
    nm706(:,i) = tmp_nm706(:,4);
    
        chrg_current(1000,i+102) = tmp_nm706(1000,3);
        chrg(:,i+102) = tmp_nm706(:,3);
        chrg_current(1:999,i+102) = 1/time_delta*diff(tmp_nm706(:,3));
    
        volt_appl(:,i+102) = tmp_nm706(:,2);
        volt_diff(1000,i+102) = tmp_nm706(1000,2);
        volt_diff(1:999,i+102) = 1/time_delta*diff(tmp_nm706(:,2));
    
    cd(loc_d);
    file_d = strcat('16Jun083',nmb,'_RTO.dat');
    disp(file_d);
    D = importdata(file_d);
    tmp_nm728 = D.data;
    
   nm728(:,i) = tmp_nm728(:,4);
   
       chrg_current(1000,i+153) = tmp_nm728(1000,3);
       chrg(:,i+153) = tmp_nm728(:,3);
       chrg_current(1:999,i+153) = 1/time_delta*diff(tmp_nm728(:,3));

        volt_appl(:,i+153) = tmp_nm728(:,2);
        volt_diff(1000,i+153) = tmp_nm728(1000,2);
        volt_diff(1:999,i+153) = 1/time_delta*diff(tmp_nm728(:,2));
end

cd(loc_e);

%% Hier das Übersichtsspetrum.

for i=1:291;
    
    nmb = num2str(i,'%03d');
    file_spectrum = strcat('16Jun09',nmb,'_RTO.dat');
    disp(file_spectrum);
    S = importdata(file_spectrum);
    tmp_spectrum = S.data;
    
    spectrum(:,i) = tmp_spectrum(:,4);

        chrg_current(1000,i+204) = tmp_spectrum(1000,3);
        chrg(:,i+204) = tmp_spectrum(:,3);
        chrg_current(1:999,i+204) = 1/time_delta*diff(tmp_spectrum(:,3));

        volt_appl(:,i+204) = tmp_spectrum(:,2);
        volt_diff(1000,i+204) = tmp_spectrum(1000,2);
        volt_diff(1:999,i+204) = 1/time_delta*diff(tmp_spectrum(:,2));
end

cd(loc_main);

chrg = chrg*c_ext;
chrg_current = c_ext*chrg_current;


%% Übersicht-Scan integrieren. Dazu einen Offset bilden zwischen 592nm und 612nm. 
    
    offset = zeros(1000);
    offset = 1/40*sum(spectrum(:,15:55),2);
    
    for i=1:291
    korrspectrum(:,i) = spectrum(:,i)-offset(:);
    end
    
    int_spectrum = zeros(291);
    int_spectrum = 1/time_delta*sum(korrspectrum(:,:),1);
    
% Daten sicher.
save('base_dat.mat');

% cd(loc_main);
% load('base_dat.mat');

%% Ersten Bilder machen. Übersichtspektrum, integriertes Spektrum und Differenzenquotienten des Entladungsstromes.

    f = figure;
    hold on;
    meshc(time_volt/1e-6,linspace(0,495,495),chrg_current'*1e3);view(2);
    xlabel('time in µs');
    ylabel('experiment #');
    % title('charge current over experiment duration');
    c = colorbar;
    c.Label.String = 'current in mA';
    box on;set(gca,'Layer','top');
%     savefig('chrg_current.fig');
    
%         saveas(gcf,'chrg_current','pdf');
        print('chrg_current','-dpdf','-noui','-bestfit');
    
    hold off;
    close(f);
    
    f=figure;
    hold on; box on;
    meshc(time_volt/1e-6,scan_spec,real(log10(korrspectrum))');view(2);
    ylabel('wavelength in nm');
    xlabel('time in µs');
    view(2);
    % title('spectrum scan via full slit');   
    c = colorbar;
    c.Label.String = 'log_{10} of intensity, a.u.';
    set(gca,'Layer','top');
    
%     savefig('spectrum.fig');
    
%         saveas(gcf,'spectrum','pdf');
        print('spectrum','-dpdf','-noui','-bestfit');
    
    hold off;
    close(f);
    
    f = figure;
    hold on; box on;
    plot(scan_spec,(-int_spectrum*1e-9+min(-int_spectrum*1e-9))/max(-int_spectrum*1e-9+min(-int_spectrum*1e-9)));
    y = -int_spectrum*1e-9/max(-int_spectrum*1e-9);
    x = scan_spec;
    axis([min(x) max(x) min(y)-0.1*max(abs(y)) max(y)+0.1*max(abs(y))]);
    xlabel('wavelength in nm');
    ylabel('norm. intensity, a.u.');
    
            % Create textbox
            annotation(f,'textbox',...
                [0.638916776595215 0.709051507397414 0.14629004648172 0.194864042981154],...
                'String',{'He I','2^3P- 3^3S','706,66 nm'},...
                'HorizontalAlignment','center',...
                'EdgeColor','none',...
                'Fontsize',18,...
                'FontName','L M Roman12');

            % Create textarrow
            annotation(f,'textarrow',[0.207907356661045 0.15008431703204],...
                [0.387284521059179 0.303625377643505],...
                'String',{'He I','2^3P-3^3D','587,65 nm'},...
                'HorizontalAlignment','center',...
                'Fontsize',18,...
                'FontName','L M Roman12');

            % Create textarrow
            annotation(f,'textarrow',[0.529761438761308 0.569983136593592],...
                [0.463696369636964 0.345921450151057],...
                'String',{'He I','2^3P-3^3S','667,98 nm'},...
                'HorizontalAlignment','center',...
                'Fontsize',18,...
                'FontName','L M Roman12');

            % Create textarrow
            annotation(f,'textarrow',[0.845314256604833 0.892074198988196],...
                [0.409760480224315 0.243202416918429],...
                'String',{'He I','2^1P - 3^3S','728,31 nm'},...
                'HorizontalAlignment','center',...
                'Fontsize',18,...
                'FontName','L M Roman12');
    
    set(gca,'Layer','top');
    % title('spectrum integrated via discharge time');
    savefig('int_spectrum');

        saveas(gcf,'int_spectrum','bmp');
%         print('int_spectrum','-dpdf','-noui','-bestfit');
    
    hold off;
    close(f);
    
    f = figure;
    hold on; box on;
    plot(scan_spec,real(log10(-int_spectrum*1e-9/max(-int_spectrum*1e-9))));
    y = real(log10(-int_spectrum*1e-9/max(-int_spectrum*1e-9)));
    x = scan_spec;
    axis([min(x) max(x) min(y)-0.1*max(abs(y)) max(y)+0.1*max(abs(y))]);
    xlabel('wavelength in nm');
    ylabel('log_{10} of norm. intensity');
    set(gca,'Layer','top');
    % title('spectrum integrated via discharge time');
%     savefig('int_spectrum');

        saveas(gcf,'log10int_spectrum','bmp');
%         print('log10int_spectrum','-dpdf','-noui','-bestfit');
    
    hold off;
    close(f);
    
    f = figure;
    hold on; box on;
    plot(volt_appl(:,1),chrg(:,1));
    x = volt_appl(:,1);
    y = chrg(:,1);
    axis([min(x) max(x) min(y)-0.1*max(abs(y)) max(y)+0.1*max(abs(y))]);
    xlabel('U_{appl}/V');
    ylabel('Q_{ext}/C');
    set(gca,'Layer','top');
    % title('applied voltage over total charge');
%     savefig('lissajous.fig');

        saveas(gcf,'lissajous','bmp');
%         print('lissajous','-dpdf','-noui','-bestfit');
    
    hold off;
    close(f);
    
%% Fit des C_tot.

    xdata = volt_appl(1:150,1);
    ydata = chrg(1:150,1);
    
    P = polyfit(xdata,ydata,1);
    c_tot = P(1);
   
% Kapazität aus der Parallelschaltung des Volumens ohne Entladung.
    
    c_par = c_tot-c_bd;
    
    volt_gap = volt_appl*(1-c_par/c_diel)-chrg./c_diel;
    current_dis = (1+c_gap/c_diel)*(chrg_current-c_tot*volt_diff);
    
    f = figure;
    hold on; box on;
    xlabel('time in µs');
    yyaxis left
    plot(time_volt/1e-6,volt_gap(:,1),'k',time_volt/1e-6,volt_appl(:,1),'k-.');
    x = time_volt/1e-6;
    ylabel('voltage in V');
    axis([min(x) max(x) -250 1250]);
    
    yyaxis right
    plot(time_volt/1e-6,current_dis(:,1)*1000,'r-');
    ylabel('current in mA');
    % title('current/appl. & gap voltage via time');
    legend('U_{gap}','U_{app}','I_{dis}');
    x = time_volt/1e-6;
    axis([min(x) max(x) -5 25]);
    zline = refline(0,0);
    zline.Color = 'k';
    zline.LineStyle = ':';
    set(gca,'Layer','top');
              
        saveas(gca,'current_dis','bmp');
%         print('current_dis','-dpdf','-noui','-bestfit');

%     savefig('current_dis.fig');

    hold off;
    close(f);
    
%% Daten der Einzelmessungen korrigieren und darstellen.

    for i=1:51;
        korr587(:,i)=nm587(:,i)-offset(:);
        korr667(:,i)=nm667(:,i)-offset(:);
        korr706(:,i)=nm706(:,i)-offset(:);
        korr728(:,i)=nm728(:,i)-offset(:);
    end

%% Bilder der Wellenlängen-Scans.

wvlgnth = [587 667 706 728];

%587nm Scan.
    i=1;
    k = wvlgnth(i);
    k = num2str(k);
    file_name = strcat('scan_',k,'.fig');
    
    f = figure;
    hold on; box on;
    meshc(time_volt/1e-6,scan1,real(log10(korr587))');
    ylabel('wavelength in nm');
    xlabel('time in µs');
    view(2);
    set(gca,'Layer','top');
    % title('scan around 587nm');   
    c = colorbar;
    c.Label.String = 'log_{10} of intensity, a.u.';
%     savefig(file_name);
    
        two = num2str(2);
        file_name = strcat('scan_',k);
%         saveas(gcf,file_name,'pdf');
%         file_name = strcat('scan_',k,two);
        print(file_name,'-dpdf','-noui','-bestfit');
    
    hold off; close(f);

    
%667nm Scan.
    i=2;
    k = wvlgnth(i);
    k = num2str(k);
    file_name = strcat('scan_',k,'.fig');
    
    f = figure;
    hold on; box on;
    meshc(time_volt/1e-6,scan2,real(log10(korr667))');
    ylabel('wavelength in nm');
    xlabel('time in µs');
    view(2);
    set(gca,'Layer','top');
    % title('scan around 667nm');   
    c = colorbar;
    c.Label.String = 'log_{10} of intensity, a.u.';
%     savefig(file_name);
    
    
        file_name = strcat('scan_',k);
%         saveas(gcf,file_name,'pdf');
%         file_name = strcat('scan_',k,two);
        print(file_name,'-dpdf','-noui','-bestfit');

    hold off; close(f);
    

%706nm Scan.
    i=3;
    k = wvlgnth(i);
    k = num2str(k);
    file_name = strcat('scan_',k,'.fig');
    
    f = figure;
    hold on; box on;
    meshc(time_volt/1e-6,scan3,real(log10(korr706))');
    ylabel('wavelength in nm');
    xlabel('time in µs');
    view(2);
    set(gca,'Layer','top');
    % title('scan around 706nm');   
    c = colorbar;
    c.Label.String = 'log_{10} of intensity, a.u.';
%     savefig(file_name);
    
        
        file_name = strcat('scan_',k);
%         saveas(gcf,file_name,'pdf');
%         file_name = strcat('scan_',k,two);
        print(file_name,'-dpdf','-noui','-bestfit');

    hold off; close(f);


%728nm Scan.
    i=4;
    k = wvlgnth(i);
    k = num2str(k);
    file_name = strcat('scan_',k,'.fig');
    
    f = figure;
    hold on; box on;
    meshc(time_volt/1e-6,scan4,real(log10(korr728))');
    ylabel('wavelength in nm');
    xlabel('time in µs');
    view(2);
    set(gca,'Layer','top');
    % title('scan around 728nm');   
    c = colorbar;
    c.Label.String = 'log_{10} of intensity, a.u.';
%     savefig(file_name);
    
        file_name = strcat('scan_',k);
%         saveas(gcf,file_name,'pdf');
%         file_name = strcat('scan_',k,two);
        print(file_name,'-dpdf','-noui','-bestfit');

    hold off; close(f);

save('2016-06-13.mat');
    
end