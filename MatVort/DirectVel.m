function V = DirectVel(diff, omega)

DEL2 = 0.25;

nrm = sqrt(DEL2 + diff(1)*diff(1) + diff(2)*diff(2) + diff(3)*diff(3));

mult = 1 ./ (4 * pi .* nrm .* nrm .* nrm);

V = mult * [omega(2)*diff(3) - omega(3)*diff(2),...
            omega(3)*diff(1) - omega(1)*diff(3),...
            omega(1)*diff(2) - omega(2)*diff(1)];
