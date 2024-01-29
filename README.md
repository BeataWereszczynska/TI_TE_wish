# TI_TE_wish
Calculating theoretical T1-weighted MRI images with arbitrary value of inversion time (TI) based on acquired Agilent MRI spin echo multi-slice experiment with inversion recovery (SEMS-IR) .fid data.

Calculating theoretical T2-weighted MRI images with arbitrary value of echo time (TE) based on acquired Agilent MRI multi-echo multi-slice (MEMS) .fid data.

The script can also be used for calculating T1, T2 and Mo maps without calculating theoretical images.


## The repository contains:
1. Python script **TI_TE_wish.py**.
2. Sample FID data in **sems_20190407_07.fid** and **mems_20190519_01.fid** folders.
3. Sample results in **Theoretical_MRI_T1w** and **Theoretical_MRI_T2w** folders.


Documentation coming soon.


## Reference: Collection of the sample data
Beata Wereszczyńska, ***Alcohol-fixed specimens for high-contrast post-mortem MRI***, Forensic Imaging, Volume 25, 2021, 200449, ISSN 2666-2256, https://doi.org/10.1016/j.fri.2021.200449. (https://www.sciencedirect.com/science/article/pii/S2666225621000208)

## License
The software is licensed under the **MIT license**. The non-software content of this project is licensed under the **Creative Commons Attribution 4.0 International license**. See the LICENSE file for license rights and limitations.

## You may also like
**MRI_k-space-derived_details_edges** - k-space based details/edges detection in MRI images with optional k-space based denoising and detail control
(data import suitable for Agilent FID files, three binarization methods to choose from), https://doi.org/10.5281/zenodo.7388435 (https://github.com/BeataWereszczynska/MRI_k-space-derived_details_edges).

**Swelling_tablet_fronts_D_k_from_MRI_T2_or_img** - Tool for characterizing the swelling of tablets immersed in a solution. Creates time plots and calculates eroding front's diffusion rate D and the rate of the swelling k from time series of either T2-maps or (properly contrasted, e.g. T2-weighted) MRI images in FDF (file format native for Agilent MRI scanners) or Text Image format. This software is suitable for swelling matrix tablets forming a membrane-like structure in contact with the solution in which they are immersed, https://doi.org/10.5281/zenodo.7262466 (https://github.com/BeataWereszczynska/Swelling_tablet_fronts_D_k_from_MRI_T2_or_img).
