function nanoscat_display_filters (psi, phi, lp)
res = 1;
figure
for j = 1:numel(psi{res})
    plot (psi{res}{j});
    hold on
end
plot (phi{res});
plot (lp, 'k')
title ('PSI/PHI  at higher resolution (and Littlewood-Paley)')
end
