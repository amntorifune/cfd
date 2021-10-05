WIP, an MUSCL version of prj1
# ns2dmuscl1
Uses staggered mesh. For the boundary cells 1st-order upwind scheme is used muscl scheme requires more cells outside the boundary to do interpolation, and I haven't yet found a way to extrapolate the values of these 'ghost cells'. For the rest of the domain, muscl scheme (3rd-order upwind) is used.
It's results are shown in figures with 'staggered half muscl' in their names.
