
This codes is for performing landmark/measure/currents/normal cycles matching in the Large Deformation Diffeomorphic Metric Mapping framework (LDDMM)

Demo programs are located in this folder :
demo_gui : graphical interactive demo (for landmark and measure matching examples only)
demo_curves_Normal_Cycles : example of matching using normal cycle metrics on curves
demo_surf_currents : example of matching using currents metrics on surfaces and Gpu implementatio for convolutions

 Notes :
 * Regarding speed, there are three possibilities for computing the kernel convolutions. This can be changed in the body of function MatchLandmarks (lines 12 to 14) : 
 - Matlab direct implementation, which is the default. This can be slower depending on the version of Matlab and/or the operating system. Typically it can be very slow on Matlab running on Linux.
 - Implementation through mex files. Mex files are stored in the MexKernels folder. When using this mode, if no correct files are found in the folder, the Matlab program will try to recompile the mex files before running.
- Implementation through Mex files that use Cuda. This is useful only for large number of points (between 10^2 and 10^5). This requires a Nvidia Gpu card with Cuda installation.
* I have included also a small demo program "script_meas" for performing unlabelled points matching using measures. Please email me if you would like to use it and want more information.

For any questions please email me :
joan.glaunes@gmail.com
 
