function stark_split3

warning off;

% Initialisierung

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


% Dateipfade der Messdaten.

loc_dat = 'E:\documents\git\ag_praktikum_2016\2016-06-28\Daten\';
loc_main = 'E:\documents\git\ag_praktikum_2016\2016-06-28\';

% loc_dat= '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-22/Daten';
% loc_main = '/Users/pha-Mac/Documents/Arbeitsgruppenpraktikum/2016-06-22';

% Beschaffe mir die Dimensionen.

% cd(loc_dat);
% tmp = importdata('16Jun60001_RTO.dat');
% tmp = tmp.data;
% [n channel] = size(tmp);
% time_volt = tmp(:,1);
% time_delta = time_volt(2)-time_volt(1);
% cd(loc_main);

% Datei-Nummern und Wellenlängen.

% m=1620;
% data_nmb = 1:1:m;
% loop_nmb = 1:1:15;
% wvl_nmb = 1:1:36;
% vert_nmb = 1:1:3;
% 
% raw_data = zeros(n,m);
% 
% rogowski = zeros(n,m);
% 
% wavelength = 491.8:0.02:492.5;
% vertical = 6.8:0.15:7.1;

% Entnehmen der Datein.
% Es kommen zuerst die 10 loops zu einer Wellenlänge. Dann erneut 10 loops
% für die nächste. Macht 360 Datein, bis eine Höhe abgearbeitet ist.

% cd(loc_dat);
% 
% for h=vert_nmb-1;
%     
%     for w=wvl_nmb-1;
%         
%         for l=loop_nmb;
% 
%             file_nmb = num2str(60000+l+15*w+540*h);
%             nmb = l+15*w+540*h;
%             disp(file_nmb);
%             file = strcat('16Jun',file_nmb,'_RTO.dat');
%             tmp = importdata(file);
%             tmp = tmp.data;
% 
%             raw_data(:,nmb) = tmp(:,channel-1);
%             
%             rogowski(:,nmb) = tmp(:,channel);
%         
%         end
%         
%     end
%     
% end

cd(loc_main);
% save('base_dat.mat');
load('base_dat.mat');

% Für refline in 3D plots.

mintime = min(time_volt*1e6);
maxtime = max(time_volt*1e6);

% Kontrolle des Experimentes
    
    f = figure;
    hold on;
    meshc(time_volt/1e-6,linspace(1,1620,1620),rogowski');view(2);
    xlabel('time in µs');
    ylabel('vertical slit pos. in inch');
    c = colorbar;
    c.Label.String = 'voltage in V';
    % title('rowgowski coil via experiment duration');
    set(gca,'YTick',[270 810 1350]);
    set(gca,'YTickLabel',{'6.8' '6.95' '7.1'});
    ax = gca;
    ax.YTickLabelRotation = 90;
    box on;set(gca,'Layer','top');
    
%     saveas(gcf,'rowgowski_full','pdf');
    print('rowgowski_full','-dpdf','-noui','-bestfit');
    
%     savefig('rowgowski_full');
    hold off; close(f);
    
f = figure;hold on;
meshc(time_volt/1e-6,wavelength,-raw_data(:,1081:15:end)'-min(min(raw_data(:,1081:15:end)')));view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
c = colorbar;
c.Label.String = 'intensity, a.u.'; 
% title('vertical position 7.1 inch');
box on;set(gca,'Layer','top');

%     saveas(gcf,'stark_71in_raw','pdf');
    print('stark_71inraw','-dpdf','-noui','-bestfit');
 
% savefig('stark_71inwar.fig');
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
    stark_695in(:,i) = 1/15*sum(raw_data(:,(i-1)*15+541:(i-1)*15+555),2);
    stark_71in(:,i) = 1/15*sum(raw_data(:,(i-1)*15+1081:(i-1)*15+1095),2);
    
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
meshc(time_volt/1e-6,wavelength,-sgf_68in'-min(min(-sgf_68in)));view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
c = colorbar;
c.Label.String = 'intensity, a.u.'; 
box on;set(gca,'Layer','top');
% title('vertical position 6.8 inch');

text(mintime + 0.25,491.85,10e-5,'6.8 inch','Color','white','FontSize',25,'FontName','L M Roman12','FontWeight','bold');

%     saveas(gcf,'stark_68in','pdf');
    print('stark_68in','-dpdf','-noui','-bestfit');
 
% savefig('stark_68in.fig');
hold off; close(f);

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,-sgf_695in'-min(min(-sgf_695in)));view(2);
ylabel('wavelength in nm');
xlabel('time in µs'); 
box on;set(gca,'Layer','top');
% title('vertical position 6.95 inch');
c = colorbar;
c.Label.String = 'intensity, a.u.';    

text(mintime + 0.25,491.85,10e-5,'6.95 inch','Color','white','FontSize',25,'FontName','L M Roman12','FontWeight','bold');

%     saveas(gcf,'stark_695in','pdf');
    print('stark_695in','-dpdf','-noui','-bestfit');

% savefig('stark_695in.fig');
hold off; close(f);

f = figure;hold on;
meshc(time_volt/1e-6,wavelength,-sgf_71in'-min(min(-sgf_71in)));view(2);
ylabel('wavelength in nm');
xlabel('time in µs');
box on;set(gca,'Layer','top');
% title('vertical position 7.1 inch');
c = colorbar;
c.Label.String = 'intensity, a.u.';   

text(mintime + 0.25,491.85,10e-5,'7.1 inch','Color','white','FontSize',25,'FontName','L M Roman12','FontWeight','bold');

%     saveas(gcf,'stark_71in','pdf');
    print('stark_71in','-dpdf','-noui','-bestfit');
    
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

    lname = strcat('lissajous',h,'.bmp');
    
    cname = strcat('currentdis',h,'.bmp');
    
    f = figure;hold on;
    set(gca,'FontSize',26); plot(volt_appl,chrg_tmp);
    xlabel('U_{appl}/V');
    ylabel('Q_{ext}/C');    
    x = volt_appl;
    y = chrg_tmp; 
    axis([min(x) max(x) min(y)-0.1*max(abs(y)) max(y)+0.1*max(abs(y))]);
    box on;set(gca,'Layer','top');
    % title('applied voltage over total charge');
    saveas(gcf,lname,'bmp');
%     print(lname,'-dpdf','-noui','-bestfit');
%     savefig(sprintf('lissajous%s.fig',h));
    hold off; close(f);
    
    f = figure;
    hold on;
    yyaxis left
    x = time_volt(1:999)/1e-6;
    set(gca,'FontSize',26); plot(time_volt(1:999)/1e-6,volt_gap(1:999),'k',time_volt(1:999)/1e-6,volt_appl(1:999),'k-.');
    xlabel('time in µs');
    ylabel('voltage in V');
    axis([min(x) max(x) -250 1250]);
    
        
    text(0.2,275,'U_{app}','FontSize',26,'FontName','L M Roman12');
    text(-0.4,1125,'U_{gap}','FontSize',26,'FontName','L M Roman12');
    text(0.3,900,'I_{dis}','Color','red','FontSize',26,'FontName','L M Roman12');
    
    yyaxis right
    set(gca,'FontSize',26); plot(time_volt(1:999)/1e-6,current_dis(1:999)*1000,'r');
    ylabel('current in mA');
    
    axis([min(x) max(x) -5 25]);
    zline = refline(0,0);
    zline.Color = 'k';
    zline.LineStyle = ':';
    
%     legend('U_{gap}','U_{app}','I_{dis}');
    box on;set(gca,'Layer','top');
    % title('current/appl. & gap voltage via time');
    saveas(gcf,cname,'bmp');
%     print(cname,'-dpdf','-noui','-bestfit');
%     savefig(sprintf('currentdis%s.fig',h));
    hold off; close(f);
    
end

% Mache die Stark-Aufspaltung.
% Hier 7.1inch.
clear trsh Ind I inmax max nmax

[trsh71, Ind71] = max(abs(stark71_korr));
inmax71 = 0;

    for i=1:36;
        
        j=Ind71(i);
        nmax71 = stark71_korr(j,i);
        
        if (abs(inmax71)<=abs(nmax71))
            nInd71 = [i j];
            inmax71 = nmax71;
        end
    end
    
    time71 = time_volt(nInd71(2))*1e6;
    
    in1max71 = 0;
    in2max71 = 0;
    starkshift71 = -stark71_korr(nInd71(2),:);
    
    for i=1:14;
        
        tmp171 = starkshift71(1+i);
        tmp271 = starkshift71(15+i);
        
        if (abs(in1max71)<=max(tmp171))
            wInd171=i+1;
            in1max71 = tmp171;
        end
        
        if (abs(in2max71)<=max(tmp271))
            wInd271=i+15;
            in2max71 = tmp271;
        end
        
    end
            
    sep71 = wavelength(wInd271)-wavelength(wInd171);
    
    fieldstrengthsq71 = -58.557+18.116*sep71+3130.96*(sep71)^2+815.6*(sep71)^3;
    fieldstrength71 = sqrt(fieldstrengthsq71);
    
    % Hier bei 6.95inch
    
    [trsh695, Ind695] = max(abs(stark695_korr));
    inmax695 = 0;

    for i=1:36;
        
        j=Ind695(i);
        nmax695 = stark695_korr(j,i);
        
        if (abs(inmax695)<=abs(nmax695))
            nInd695 = [i j];
            inmax695 = nmax695;
        end
    end
    
    time695 = time_volt(nInd695(2))*1e6;
    
    in1max695 = 0;
    in2max695 = 0;
    starkshift695 = -stark695_korr(nInd695(2),:);
    
    for i=1:14;
        
        tmp1695 = starkshift695(1+i);
        tmp2695 = starkshift695(15+i);
        
        if (abs(in1max695)<=max(tmp1695))
            wInd1695=i+1;
            in1max695 = tmp1695;
        end
        
        if (abs(in2max695)<=max(tmp2695))
            wInd2695=i+15;
            in2max695 = tmp2695;
        end
        
    end
            
    sep695 = wavelength(wInd2695)-wavelength(wInd1695);
    
    fieldstrengthsq695 = -58.557+18.116*sep695+3130.96*(sep695)^2+815.6*(sep695)^3;
    fieldstrength695 = sqrt(fieldstrengthsq695);
    
    % Hier bei 6.8inch
    
    [trsh68, Ind68] = max(abs(stark68_korr));
    inmax68 = 0;

    for i=1:36;
        
        j=Ind68(i);
        nmax68 = stark68_korr(j,i);
        
        if (abs(inmax68)<=abs(nmax68))
            nInd68 = [i j];
            inmax68 = nmax68;
        end
    end
    
    time68 = time_volt(nInd68(2))*1e6;
    
    in1max68 = 0;
    in2max68 = 0;
    starkshift68 = -stark68_korr(nInd68(2),:);
    
    for i=1:14;
        
        tmp168 = starkshift68(1+i);
        tmp268 = starkshift68(15+i);
        
        if (abs(in1max68)<=max(tmp168))
            wInd168=i+1;
            in1max68 = tmp168;
        end
        
        if (abs(in2max68)<=max(tmp268))
            wInd268=i+15;
            in2max68 = tmp268;
        end
        
    end
            
    sep68 = wavelength(wInd268)-wavelength(wInd168);
    
    fieldstrengthsq68 = -58.557+18.116*sep68+3130.96*(sep68)^2+815.6*(sep68)^3;
    fieldstrength68 = sqrt(fieldstrengthsq68);
    
    f = figure;hold on;
    set(gca,'FontSize',26); plot(wavelength,starkshift68-min(starkshift68),wavelength,starkshift695-min(starkshift695),wavelength,starkshift71-min(starkshift71));
    x = wavelength;
    ymax = starkshift71-min(starkshift71);
    ymin = starkshift68-min(starkshift68);
    axis([min(x) max(x) min(ymin)-0.1*max(abs(ymin)) max(ymax)+0.1*max(abs(ymax))]);
    xlabel('wavelength in nm');
    ylabel('intensity, a.u.');
    
        line1 = line([wavelength(wInd171) wavelength(wInd171)],[min(ymin)-0.1*max(abs(ymin)) max(ymax)+0.1*max(abs(ymax))]);
        line1.Color = 'black';
        line1.LineStyle = ':';
    
        line2 = line([wavelength(wInd271) wavelength(wInd271)],[min(ymin)-0.1*max(abs(ymin)) max(ymax)+0.1*max(abs(ymax))]);
        line2.Color = 'black';
        line2.LineStyle = ':';
        
        tmp = num2str(wavelength(wInd171)); wvl171 = strcat(tmp,' nm'); text(wavelength(wInd171)+0.01,17e-5,wvl171,'Color','black','FontSize',22,'FontName','L M Roman12');
        tmp = num2str(wavelength(wInd271)); wvl271 = strcat(tmp,' nm'); text(wavelength(wInd271)+0.01,17e-5,wvl271,'Color','black','FontSize',22,'FontName','L M Roman12');
       
        annotation(f,'doublearrow',[0.55859375 0.365234375],[0.84765625 0.84765625]);
        text('Position',[492.032115869018 0.000151581873091233 0],'String','\Delta\lambda = 0.2 nm','Color','black','FontSize',22,'FontName','L M Roman12');
        
        text('Position',[492.27 11.5e-5],'String','7,1 inch','Color','blue','FontSize',22,'FontName','L M Roman12');
        text('Position',[492.01 3e-5],'String','6,95 inch','Color','red','FontSize',22,'FontName','L M Roman12');
        text('Position',[492.1 0.5e-5],'String','6,8 inch','Color','black','FontSize',22,'FontName','L M Roman12');
        
    box on;set(gca,'Layer','top');
    % title(sprintf('shift, 7.1in, %gµs, E=%g kV/cm', time, fieldstrength));
    saveas(gcf,'stark_shift','bmp');
%     print('stark_shift71in','-dpdf','-noui','-bestfit');
%     savefig(sprintf('stark_shift71in',h));
    hold off; close(f);
    
%% Daten.

save('2016-06-28.mat');

end