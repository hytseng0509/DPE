## Direct 3D Pose Estimation of a Planar Target (WACV 2016)

### Introduction

This is the research code for the paper:

Hung-Yu Tseng, Po-Chen Wu, Ming-Hsuan Yang, and Shao-Yi Chien, "Direct 3D Pose Estimation of a Planar Target", WACV 2016

The proposed algorithm estimates the 3D pose of a planar target with respect to a calibratede camera.
It performs favorably against state-of-the art feature-based methods in terms of accuracy and robustness.

### Citation

Please cite the paper if you find the code useful in your research.

    @inproceedings{tsengdirect,
        title={Direct 3D Pose Estimation of a Planar Target},
        author={Tseng, Hung-Yu and Wu, Po-Chen and Yang, Ming-Hsuan and Chien, Shao-Yi},
        booktitle = {Proc. IEEE Winter Conference on Applications of Computer Vision},
        year = {2016}
    }
    
### Contents
|Name|Descriptions|
|---|---|
|imgs|camera image and target image for the demo|
|matlab|matlab codes|
|mex|mex function codes|
|Test_DPE.m|for running the whole proposed algorithm|
|Test_APE.m|for running the proposed approximated pose estimation|
|Test_Refine.m|for running the proposed refinement scheme|

To run the whole algorithm, run the "Demo.m".
To run the approximated pose estimation, comment line 26 and un-comment line 27 in the "Demo.m" and run it.
To run the proposed refinement scheme, you need to provide an initial extrinsic matrix.

Note that some of the parameters have been optimized compared to the version used in the paper.
The performance may be slightly different between these two versions.

If you have any problem, please contact Hung-Yu Tseng.

### Reference
Y. Zheng, Y. Kuang, S. Sugimoto, K. Astrom, and M. Okutomi, “Revisiting the PnP Problem: A Fast, General and Optimal Solution,” in Proc. IEEE International Conference on Computer Vision, 2013.