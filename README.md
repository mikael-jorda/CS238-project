# Real Time Collision Detection and Identification for Robotic Manipulators

Real time collision detection and identification on a simulated 4-DOF planar robot. In this repo you will find the implementation of a:
  - **particle filter** - probabilistic method
  - **momentum observer** - analytical method
  
The performance of these two methods in terms of force estimation and localization is studied during the excution of two tasks: 
  - holding position 
  - tracking a trajectory

Below is a figure depicting the 4-DOF robot and the PF convergence. The green line is the track trajectory and the blue sphere the object that is colliding with the robot. Further details and experimental results can be found here: [**ArXiv Paper**](https://arxiv.org/abs/1802.00546)

<br />
<p align="center">
  <img src="https://github.com/mikael-jorda/CS238-project/blob/master/readme_fig.PNG" width="800"/>
</p>
<br />
