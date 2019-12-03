function [h_colorbar,h_surf,h_contour,fig1] = surf_and_con(XI,YI,ZI,fig1)
    isFigureHandle = 0;
    
    if (exist('fig1','var'))
%         isFigureHandle = ishandle(fig1) & strcmp(get(fig1,'type'),'figure');
       isFigureHandle = any(ishandle(fig1));
    end
    
    if (~exist('fig1','var') || isFigureHandle == 0)
        figure1= figure('PaperType','A4','Visible','on','PaperPositionMode', 'auto');
        fig1 = axes('Parent',figure1,'Layer','top','FontSize',16);
        hold all
        grid on
        box on
    end

    h_surf = surf(fig1,XI,YI,ZI,'Linestyle','none');
    
    h_colorbar = colorbar('peer',fig1);

    set(get(h_colorbar,'ylabel'),'Fontsize', 16);
    set(h_colorbar,'Fontsize',16);

    if (max(max(ZI)) > 0 && min(min(ZI)) >= 0)
        set(fig1,'ZDir','reverse');
    end
    
    if (max(max(ZI)) > 0 && min(min(ZI)) < 0)
        delete(h_surf);
        delete(h_colorbar);
        scale_factor = 10^(-floor(log10(abs(min(min(ZI)))))+1);
        h_surf = surf(fig1,XI,YI,ZI-scale_factor,'Linestyle','none');
        h_colorbar = colorbar('peer',fig1);
        set(get(h_colorbar,'ylabel'),'Fontsize', 16);
        set(h_colorbar,'Fontsize',16);
        ytick_lab = get(h_colorbar,'yticklabel');
        for j = 1:length(ytick_lab)
            ytick_lab_new{j} = str2num(ytick_lab{j})+scale_factor;
        end
        set(h_colorbar,'yticklabel',ytick_lab_new);
        caxis_range = caxis(fig1);
    end

    shading(fig1,'interp');

    xlim(fig1,[min(min(XI)) max(max(XI))]);
    ylim(fig1,[min(min(YI)) max(max(YI))]);
    view(fig1,[0,90]);
    
    [C,h_contour]= contour(fig1,XI,YI,ZI,'Color','k');
    clabel(C,h_contour); 
%     clabel(C,h,'Fontsize',14); 
    
    if (max(max(ZI)) > 0 && min(min(ZI)) < 0)
        caxis(fig1,caxis_range);
    end
end