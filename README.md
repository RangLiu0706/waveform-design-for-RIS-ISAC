### About the paper
This is a code package for the paper: 
R. Liu, M. Li, Y. Liu, Q. Wu, and Q. Liu, “Joint transmit waveform and passive beamforming design for RIS-aided DFRC systems,”IEEE J. Sel. Topics Signal Process., vol. 16, no. 5, pp. 995-1010, Aug. 2022.

@ARTICLE{9769997,
  author={Liu, Rang and Li, Ming and Liu, Yang and Wu, Qingqing and Liu, Qian},
  journal={IEEE Journal of Selected Topics in Signal Processing}, 
  title={Joint Transmit Waveform and Passive Beamforming Design for RIS-Aided DFRC Systems}, 
  year={2022},
  volume={16},
  number={5},
  pages={995-1010},
  keywords={Radar;Sensors;Clutter;Wireless communication;Quality of service;Array signal processing;Signal to noise ratio;Reconfigurable intelligent surface (RIS);dual-functional radar-communication (DFRC);multi-input multi-output (MIMO) radar;multi-user multi-input single-output (MU-MISO);space-time adaptive processing (STAP)},
  doi={10.1109/JSTSP.2022.3172788}}


- If you use this simulation code package in any way, please cite the original paper above.
- All codes are contributed by Rang Liu (email: rangl2@uci.edu; website: https://rangliu0706.github.io/). 
   Please feel free to contact with her if you have any suggestions. 
- The link of this paper is: https://ieeexplore.ieee.org/document/9769997
- More information can be found at: https://www.minglabdut.com/resource.html
- Copyright Notice: This code is licensed for personal, non-commercial use only, specifically for academic purposes. Copyright reserved by the MingLab (led by Prof. Ming Li), School of Information and Communication Engineering, Dalian University of Technology, Dalian 116024, China. 


### Software platform
- Please note that the MATLAB2022b is used for this simulation code package, and there may be some imcompatibility problems among different sofrware versions. 
- To run those codes, please download and install [CVX](http://cvxr.com/cvx/) & [Manopt](https://www.manopt.org/)

### Content of this simulation code package
- The files "main_iterations", "main_SINR_P", "main_SINR_N", "main_SINR_SNR", "main_SINR_M", "main_N_alpha", and "main_SINR_d" are used to generate Figs. 3-9, respectively.

Abstract of the paper: 
Reconfigurable intelligent surface (RIS) is a promising technology for 6 G networks owing to its superior ability to enhance the capacity and coverage of wireless communications by smartly creating a favorable propagation environment. In this paper, we investigate the potential of employing RIS in dual-functional radar-communication (DFRC) systems for improving both radar sensing and communication functionalities. In particular, we consider a RIS-assisted DFRC system in which the multi-antenna base station (BS) simultaneously performs both multi-input multi-output (MIMO) radar sensing and multi-user multi-input single-output (MU-MISO) communications using the same hardware platform. We aim to jointly design the dual-functional transmit waveform and the passive beamforming of RIS to maximize the radar output signal-to-interference-plus-noise ratio (SINR) achieved by space-time adaptive processing (STAP), while satisfying the communication quality-of-service (QoS) requirement under one of three metrics, the constant-modulus constraint on the transmit waveform, and the unit-modulus constraint of RIS reflecting coefficients. An efficient algorithm framework based on the alternative direction method of multipliers (ADMM) and majorization-minimization (MM) methods is developed to solve the complicated non-convex optimization problem. Simulation results verify the advancement of the proposed RIS-assisted DRFC scheme and the effectiveness of the developed ADMM-MM-based joint transmit waveform and passive beamforming design algorithm.



