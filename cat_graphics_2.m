clear all
clc
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


colorType = 0; %0 for regular, 1 for F=3,4 block colors
make_movie = 1;
saveImages = 1;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%     rfxWave = mod(opt_params.rf_wave(1,:),2*pi);
%     rfyWave = mod(opt_params.rf_wave(2,:),2*pi);
%     mwWave = mod(opt_params.mw_wave(1,:),2*pi);

    

    %%AARON
%    times = linspace(0,210,wavePoints);
  
    




psi=zeros(10,1);
psi(2)=1/sqrt(2);
psi(9)=1/sqrt(2);
rho_target=psi*psi'



% opt_params.rf_wave = zeros(2,num_rhos-1);
% opt_params.rf_wave(1,:) = opt_params.control_fields_final(1:num_rhos-1,1);
% opt_params.rf_wave(2,:) = opt_params.control_fields_final(1:num_rhos-1,2);
% 
% opt_params.mw_wave = zeros(1,num_rhos-1);
% opt_params.mw_wave(:) = opt_params.control_fields_final(1:num_rhos-1,3);
% 
% 
% dt = 1e-6;
% %resample waveforms every microsecond
% opt_params.points = (num_rhos-1)*5;
% opt_params.rf_wave = [rectpulse(opt_params.rf_wave(1,:),5);rectpulse(opt_params.rf_wave(2,:),5)];
% opt_params.mw_wave = rectpulse(opt_params.mw_wave(1,:),5);
% 
% points = opt_params.points + 1;
% 
% opt_params.control_fields = zeros(opt_params.points,3);
% opt_params.control_fields(:,1) = opt_params.rf_wave(1,:);
% opt_params.control_fields(:,2) = opt_params.rf_wave(2,:);
% opt_params.control_fields(:,3) = opt_params.mw_wave(:);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Open new figure with specific size
figure('Color',[1 1 1],'Units','inches','OuterPosition',[1,2,12.5,9.5])


%%%%%
%%%%  target 1
%%%%

%% MAKE SMART COLOR BAR CHART OF TARGET MATRIX
%subplot(5,2,[5 7 9])
subplot(2,2,3)
barHandlesTarget = bar3(abs(rho_target));
%text(19,-11.5,'Target State','FontSize',22,'FontName','Times');
% text('Interpreter','latex',...
%     'String','$$\rho$$ target',...
%     'Position',[19 -11.5],'FontSize',24,'Color','k');
% text('Interpreter','latex','String',...
%     '$$|\psi_{\rm target}\rangle=|3,0\rangle$$',...
% 	'Position',[16.5 -7.5],'FontSize',14);
% text('Interpreter','latex','String',...
%     '$$|\psi_{\rm target}\rangle=\frac{1}{\sqrt{2}}\big(|3,3\rangle+|3,-3\rangle\big)$$',...
% 	'Position',[11.5 -9.5],'FontSize',14);
%text('Interpreter','latex','String',...
%    '$$|\psi_{\rm target}\rangle=\big(|4,-2\rangle\big)$$',...
%	'Position',[10.5 -11],'FontSize',14);

colorMapHandles = get(gcf,'colormap');  % Use the current colormap.
count = 0;

for k = 1:length(barHandlesTarget)
    xd = get(barHandlesTarget(k),'xdata');
    yd = get(barHandlesTarget(k),'ydata');
    zd = get(barHandlesTarget(k),'zdata');
    delete(barHandlesTarget(k))
    idx = [0;find(all(isnan(xd),2))];
    
    if k == 1
        barSurf = zeros(length(barHandlesTarget)*(length(idx)-1),1);
        dv = floor(size(colorMapHandles,1)/length(barSurf));
    end
    
    for q = 1:length(idx)-1
        count = count + 1;
        barSurf(count) = surface(xd(idx(q)+1:idx(q+1)-1,:),...
            yd(idx(q)+1:idx(q+1)-1,:),...
            zd(idx(q)+1:idx(q+1)-1,:),...
            'facecolor',colorMapHandles((count-1)*dv+1,:));
    end
end

%Define each element (bar) of the dim^2 plot (barSurf) to be part of one of
%the diagonals.  Diagonals are defined below:
% |:  :  :  :  :     |
% |                  |
% |7  5  3  1  2  ...|
% |                  |
% |5  3  1  2  4  ...|
% |                  |
% |3  1  2  4  6  ...|
% |                  |
% |1  2  4  6  8  ...|

%Figure out which elements are in each diagonal
dim = 10;
setNumbers = ones(dim*dim,1);
for k = 1:dim
    indexOne = (k-1)*dim + 1;
    indexTwo = k*dim;
    setNumbers(indexOne:indexTwo,1) = [fliplr(1:2:2*k-1),2:2:2*(dim-k)];
end

%Figure out which elements are F=3,4 or both


% Define the color for each of the diagonals
colorFactors = [1,rectpulse(linspace(1,1,dim-1),2)];
colors = zeros(2*dim-1,3);
colors(1,:) = [0.80 0.08 0.08];
colors(2,:) = [0.60 0.80 1.00];
colors(3,:) = [0.60 0.80 1.00];
colors(4,:) = [0.60 0.80 1.00];
colors(5,:) = [0.60 0.80 1.00];
colors(6,:) = [0.60 0.80 1.00];
colors(7,:) = [0.60 0.80 1.00];
colors(8,:) = [0.60 0.80 1.00];
colors(9,:) = [0.60 0.80 1.00];
colors(10,:) = [0.60 0.80 1.00];
colors(11,:) = [0.60 0.80 1.00];
colors(12,:) = [0.60 0.80 1.00];
colors(13,:) = [0.60 0.80 1.00];
colors(14,:) = [0.60 0.80 1.00];
colors(15,:) = [0.60 0.80 1.00];
colors(16,:) = [0.60 0.80 1.00];
colors(17,:) = [0.60 0.80 1.00];
colors(18,:) = [0.60 0.80 1.00];
colors(19,:) = [0.60 0.80 1.00];
colors(20,:) = [0.60 0.80 1.00];
colors(21,:) = [0.60 0.80 1.00];
colors(22,:) = [0.60 0.80 1.00];
colors(23,:) = [0.60 0.80 1.00];
colors(24,:) = [0.60 0.80 1.00];
colors(25,:) = [0.60 0.80 1.00];
colors(26,:) = [0.60 0.80 1.00];
colors(27,:) = [0.60 0.80 1.00];
colors(28,:) = [0.60 0.80 1.00];
colors(29,:) = [0.60 0.80 1.00];
colors(30,:) = [0.60 0.80 1.00];
colors(31,:) = [0.60 0.80 1.00];

%Loop through the bars and set the color of each bar
for k = 1:100
    setNum = setNumbers(k,1);
    if colorType == 0
        set(barSurf(k),'facecolor',colors(setNum,:)*colorFactors(setNum));
    elseif colorType ==1
        if max(k == f4elements) == 1
            colorVec = [0.60 0.80 1.00];
            if setNum == 1
                colorVec = [0.80 0.08 0.08];
            end
        elseif max(k == f3elements) == 1
            colorVec = [0.60 0.70 1.00];
            if setNum == 1
                colorVec = [0.80 0.3 0.08];
            end
        else
            colorVec = [0.40 0.60 0.80];
        end
        set(barSurf(k),'facecolor',colorVec);
    end
end


zlim([0 1])
currentAxes = get(gcf,'CurrentAxes');
set(currentAxes,'XGrid','off');
set(currentAxes,'YGrid','off');
set(currentAxes,'ZGrid','off');

set(currentAxes,'ZTick',[0,0.5,1]);
  set(currentAxes,'XTick',[1,10,16]);
    set(currentAxes,'YTick',[1,10,16]);
    set(currentAxes,'XTickLabel',{'|9/2,-9/2>';'|9/2,9/2>';'|3,-3>'},'FontSize',14)
    set(currentAxes,'YTickLabel',{'|9/2,-9/2>';'|9/2,9/2>';'|3,-3>'},'FontSize',14)
mwWave=importdata("D:\Graphics_state_prep\final_cat_wave.txt");
mwWave=[mwWave(1);mwWave];

B1=importdata("D:\Graphics_state_prep\final_cat_real_"+5+".txt");
B2=importdata("D:\Graphics_state_prep\final_cat_imag_"+5+".txt")
J=10
   for i=1:151
       F(i)=2*(i-1)/150;
   end
for i=1:151
   i
    %time=2*(i-1)/30;
    psi=B1(1:J,i)+1i*B2(1:J,i);
    rho=psi*ctranspose(psi);
     subplot(2,2,1:2)
    stairs(F,mwWave,'LineWidth',2);
    hold on
    plot(F(i),mwWave(i),'o','MarkerEdgeColor','r',...
         'MarkerFaceColor','b','MarkerSize',10);
     hold off
     xlabel('$\Omega T/\pi$','Interpreter','latex')
     ylabel('$C(t)/\pi$','Interpreter','latex')
     set(gca,'Fontsize',16)
    subplot(2,2,4)

      barHandles = bar3(abs(rho));
    %text(18,-11,'Evolved State','FontSize',22,'FontName','Times');
%     text('Interpreter','latex',...
%         'String','$$\rho$$ evolved',...
%         'Position',[19 -11.5],'FontSize',24,'Color','k');
    
    colorMapHandles = get(gcf,'colormap');  % Use the current colormap.
    count = 0;
    
    for k = 1:length(barHandles)
        xd = get(barHandles(k),'xdata');
        yd = get(barHandles(k),'ydata');
        zd = get(barHandles(k),'zdata');
        delete(barHandles(k))
        idx = [0;find(all(isnan(xd),2))];
        
        if k == 1
            barSurf = zeros(length(barHandles)*(length(idx)-1),1);
            dv = floor(size(colorMapHandles,1)/length(barSurf));
        end
        
        for q = 1:length(idx)-1
            count = count + 1;
            barSurf(count) = surface(xd(idx(q)+1:idx(q+1)-1,:),...
                yd(idx(q)+1:idx(q+1)-1,:),...
                zd(idx(q)+1:idx(q+1)-1,:),...
                'facecolor',colorMapHandles((count-1)*dv+1,:));
        end
    end
    
    dim = 10;
    setNumbers = ones(dim*dim,1);
    for k = 1:dim
        indexOne = (k-1)*dim + 1;
        indexTwo = k*dim;
        setNumbers(indexOne:indexTwo,1) = [fliplr(1:2:2*k-1),2:2:2*(dim-k)];
    end
    
    %Loop through the bars and set the color of each bar
    for k = 1:100
        setNum = setNumbers(k,1);
        if colorType == 0
            set(barSurf(k),'facecolor',colors(setNum,:)*colorFactors(setNum));
        elseif colorType ==1
            if max(k == f4elements) == 1
                colorVec = [0.60 0.80 1.00];
                if setNum == 1
                    colorVec = [0.80 0.08 0.08];
                end
            elseif max(k == f3elements) == 1
                colorVec = [0.60 0.70 1.00];
                if setNum == 1
                    colorVec = [0.80 0.3 0.08];
                end
            else
                colorVec = [0.40 0.60 0.80];
            end
            set(barSurf(k),'facecolor',colorVec);
        end
    end
    
   
    zlim([0 1])
    currentAxes = get(gcf,'CurrentAxes');
    set(currentAxes,'XGrid','off');
    set(currentAxes,'YGrid','off');
    set(currentAxes,'ZGrid','off');
   
%      subplot(8,2,[1 2 3 4])
%      plot(mwWave,'o','MarkerEdgeColor','r',...
%         'MarkerFaceColor','b','MarkerSize',10);
  
 
 
    zlim([0,1])
   
    set(currentAxes,'XTick',[1,10,16]);
    set(currentAxes,'YTick',[1,10,16]);
    set(currentAxes,'XTickLabel',{'|9/2,-9/2>';'|9/2,9/2>';'|3,-3>'},'FontSize',14)
    set(currentAxes,'YTickLabel',{'|9/2,-9/2>';'|9/2,9/2>';'|3,-3>'},'FontSize',14)
     movieVector(i)=getframe(gcf);
    
end

myWriter=VideoWriter("cat_state_graphics");
myWriter.FrameRate=10;

open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);



    
  
