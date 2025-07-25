# StabVSP
Python script for evaluate stablity using OpenVSP's stab file. Current version: v2.

# Scripts' Description
There are two scripts and can be launched separately.

The first script `Inertia.py` calculates the moment of inertias. It reads input file `Mass_Inertia_imp.txt` and outputs the total moment of inertia into `Mass_Inertia_out.txt`. Each line includes `CGx CGy CGz m(mass) Ixx Iyy Izz Ixz` in SI units. Ixy and Iyz are not calculated because the aircraft is considered to be symmetrical in these two directions when analysing the stability. 

The second script `StabVSP.ipynb` reads `.stab` file and inertia file in the specified dirs (corresponding with each configuration), saves these data in a python class, and evaluates both static and dynamic stability of the aircraft. 

For static stability, the inertia file is not needed. Just comment lines for reading inertia, and plot the aerodynamic coefficients or derivatives, such as `Cm` or `Cnbeta`.

For dynamic stability, the eigenvalues and eigenvectors of small pertubration equations in **body** coordinate system are calculated, in both lonitudinal and lateral directions. Note that this program needs the aircraft's pitch angle $\theta$ (the angle between x axis in body axes and the ground) and climbing angle $\gamma$ (the angle between velocity and the ground) as inputs. the pitch angle is used in lateral dynamic stability's calculation, and the climbing angle is used in longitudinal dynamic stability's calculation. You can set them as zero in the initial analysis.
Analysis in **wind** coordinates can also work, but the coefficient matrix needs to be rewritten. 

In version 2, the dynamic derivatives in unsteady scenes (p, q, and r analysis) will be considered. The script can read `.pstab`, `.qstab` and `.rstab` file (if exist) and fetch derivatives with `alpha_dot` and `beta_dot`. `Cm_alpha_dot` & `CL_alpha_dot` is used for longitudinal analysis, and other derivatives are not used currently. You can switch to a more approximated way (ignore `CL_alpha_dot` & related `Z_alpha_dot`) by uncomment & comment different expression of matrix `A1`.
Tips: Check the convergence and add wake iteration numbers for precision.


# Coordinate Systems

It's worth noting that the coordinate systems' direction is very important. There are four coordinate systems, and the coordinate systems used for stability analysis are according to ISO standards.

1. Ground Coordinates.
2. Aircraft Body Coordinates (Used in stability analysis): X - forward, Y - right, Z - downward. 
3. Aircraft Body Coordinates In OpenVSP: X - backward, Y - right, Z - upward. 
4. Wind Coordinates: Alpha positive if velocity vector points downwards, Beta positive if velocity vector points right.

Values in `.stab` file (CL, CD, CS, CMl, CMm, CMn, and those derivatives with alpha, beta, p, q, r) are according to coordinates 2 & 4 (ISO standard) while original data (CFx, CFz, CMx, CMy, CMz and those derivatives with alpha, beta, p, q, r) are according to coordinates 3 & 4 and not used for calculation in the python script except those related with `CFy`. That does not matter since y-direction is the same between coordinates 2 & 3.

# Usage, Examples and Standard Models

Method 1 ( v1 & v2 ) : Use `StabVSP_v*.ipynb`.

Firstly, calculate aerodynamic derivatives in OpenVSP, and inertia using `Inertia.py`(optional if you only evaluate static stability). Then drop the `.stab` file and inertia file into one directory, and write the path into the variables in `StabVSP.ipynb`. Then run the `.ipynb` script. Skip the parts about inertia files if you only evaluate static stability.

See `Examples` directory which includes a flying wing aircraft `Progress 6`. Please refer to `.vsp3` and `.vspaero` file for the settings for stablity analysis settings in OpenVSP.

Method 2 ( v2 ): copy `stabvsp.py` and use it as a package instead of writing all the things in to the `.ipynb` file. See `demo.ipynb` for usage example.

A standard model in NASA report TM-4640 has been added to validate the accuracy of VSPAERO and the program. The detailed comparison result of the aerodynamic derivatives between OpenVSP's VLM methods, VLM513's VLM methods and windtunnel test results are needed to be added. OpenVSP 3.43.0's result is not bad roughly seeing, but needed to be inspected. OpenVSP 3.45.x's result may differ since the solver has been modified a lot.
[VLM513](https://shi.buaa.edu.cn/songlei/zh_CN/jxzy/20673/content/1167.htm) is a MATLAB program developed by SONG Lei in BUAA which is more precise in calculating lateral aerodynamic derivatives. 

# Reference & Cite

This script is used in the paper “李志锴,魏莎,丁虎,等.平直翼飞翼布局飞机的操稳特性[J].上海大学学报(自然科学版),2024,30(05):925-937. ”

The small pertubration equations, the meaning of each `aerodynamic derivatives` and `state derivatives` can be found in "方振平, 陈万春, 张曙光 编著. 航空飞行器飞行动力学[M]. 北京: 北京航空航天大学出版社, 2005." The calculation of moment of inertia is based on "郭茂政. 论惯性积的平移变换和旋转变换[J]. 大学物理, 2004, 23(6): 23."

If you need to cite this program, please cite this GitHub page or cite the paper above. 
```
@misc{StabVSP,
  author = {LI Zhikai},
  title = {{StabVSP}: Python script for evaluate stablity using VSPAero's aerodynamic derivatives result. },
  year = {2025},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/Kai2510/StabVSP}},
}
