# growth

This repo contains the growth codes (Abaqus UMAT subroutines) described in the [Hitchhiker's Guide to Abaqus](https://zenodo.org/record/1243270).  If you find these helpful, please include the Hitchhiker's Guide in your attributions.

Simple, geometrically-relevant input files are included for each UMAT, such that the pairs listed in the table below should run together.  

| Input files                               | UMAT files           |
|-------------------------------------------|----------------------|
| cube_1_C3D8_noload.inp                    | umat_fiber_morph.f   |
| cube_1_C3D8_stretch_y.inp                 | umat_fiber_stretch.f |
| cube_1_C3D8_stretch_z.inp                 | umat_fiber_stretch.f |
| cube_1_C3D20_stretch_x.inp                | umat_fiber_stretch.f |
| cube_1_C3D20R_stretch_x.inp               | umat_fiber_stretch.f |
| cube_8_C3D8_stretch_x.inp                 | umat_fiber_stretch.f |
| cube_8_C3D20R_stretch_x.inp               | umat_fiber_stretch.f |
| muscle_stretch_z.inp                      | umat_fiber_stretch.f |
| square_noload.inp                         | umat_area_morph.f    |
| sheet_noload.inp                          | umat_area_morph.f    |
| circle_pressure.inp                       | umat_area_stretch.f  |
| cube_1_C3D8_stretch_xy.inp                | umat_area_stretch.f  |
| cube_8_C3D8_stretch_xy.inp                | umat_area_stretch.f  |
| muffin_noload.inp                         | umat_iso_morph.f     |
| cube_1_C3D8_stretch_xyz_iso_stretch.inp   | umat_iso_stretch.f   |
| cube_1_C3D8_stretch_xyz_iso_Mandel.inp    | umat_iso_Mandel.f    |
| cube_1_C3D8_stretch_xyz_transverse.inp    | umat_transverse.f    |
| cube_1_C3D8_stretch_xyz_ortho_stretch.inp | umat_ortho_stretch.f |

[![DOI](https://zenodo.org/badge/146362714.svg)](https://zenodo.org/badge/latestdoi/146362714)
