function [figure1,axes1] = getFigure( cascade, w, h )

if ~exist('w','var')
    w = 560;
    h = 520;
end

figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
box(axes1,'on');
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16,'FontWeight','bold','LineWidth',2);

if exist('cascade', 'var' ) 
    if cascade
        figN = get(gcf,'Number');
        
        r = groot;

        w_max = r.ScreenSize(3);
        h_max = r.ScreenSize(4);
        %w_max = 3840;
        %h_max = 2160;

        n_W = floor( w_max / w );
        n_H = floor( h_max / h );

        nMax = n_H * n_W;

        %find the position of the current window within the allowed max number
        %of windows on the screen. If figN is greater than  nMax then we circle
        %back
        figN = mod( figN, nMax );
        if figN == 0
            figN = nMax;
        end
        
        j = ceil( figN / n_W );
        i = figN - (j-1)*n_W;
        
        figure1.OuterPosition = [(i-1)*w,h_max-j*h,w,h];
    end
    
    
end

end