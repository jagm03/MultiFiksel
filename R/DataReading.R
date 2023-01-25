#Preparing data
PP <- rthin(unmark(Cells$pp$"13") %mark% Cells$pp$'13'$marks$Immune, .1)
plot(PP)

Gamma <- matrix(rep(2, 25), 5, 5)
Radii <- matrix(rep(12, 25), 5, 5)

A <- ppm(PP ~ 1, MultiFiksel(iradii = Radii, igammaii = Gamma))
plot(fitin(A))
B <- rmhmodel() 
