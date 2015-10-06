function nanoscat_display_filters (psi, phi, lp)
res = 1;
hold on
for j = 1:numel(psi{res})
    plot (psi{res}{j});
end
plot (phi{res});
plot (lp, 'k')
title ('PSI/PHI  at higher resolution (and Littlewood-Paley)');
hold off
end
