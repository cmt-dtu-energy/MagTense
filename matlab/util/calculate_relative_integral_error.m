function rel_int_error = calculate_relative_integral_error(x_f,f,x_g,g)
    
    % Interpolate g(x_g) to x_f calculate the relative integral in percent between the two curves 
    
    indx = ~isnan(g);
    g_interpolated(:,1) = interp1(x_g(indx),g(indx),x_f,'linear','extrap');
    
    rel_int_error(1) = trapz(x_f,abs(f-g_interpolated))./trapz(x_f,abs(f))*100;
end