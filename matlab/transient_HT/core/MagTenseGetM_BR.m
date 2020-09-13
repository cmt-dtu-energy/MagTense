
%By Kaspar Kirstein Nielsen, kasparkn@gmail.com, June 2020
%Finds the magnetization as a function of field, H, and temperature, T, and
%pressure field, p, given the hysteresis field, hyst (binary) and the state function, stateFct
%the input ind holds the indices into the state function for each tile
%(primitive) in the model, i.e. length(ind)=length(T).
%given in some version of the mean field / Bean Rodbell approximation
function M = MagTenseGetM_BR( H, T, p, hyst, ind, stateFct )

%we assume the magnet to be soft in the sense that M = M0 * H/Hnorm, i.e.
%the magnetization is parallel to H and has the magnitude M0. This
%magnitude we can find from the Bean-Rodbell / mean field model.

%norm of the field
Hnorm = sqrt(sum(H.^2,2));

Mnorm = stateFct.getMnorm( Hnorm, T, p, hyst, ind );



M = Mnorm .* H./Hnorm;

end