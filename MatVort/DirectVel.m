function V = DirectVel(diff, omega)

DEL = 0.25;



nrm = sqrt(DEL + diff.*diff);
mult = 1 ./ (4 * pi .* nrm .* nrm .* nrm);



V = mult.*diff.*omega;
