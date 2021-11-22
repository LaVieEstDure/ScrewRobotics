# ScrewRobotics
---
A Robotics library for control, planning and estimation algorithms that use Screw and Lie theory. Made for UQCS Hackathon 2021.

![screwgif](/_static/screw.gif)


---
# How to use
## Installation
```bash
git clone https://github.com/LaVieEstDure/ScrewRobotics/
cd ScrewRobotics
pip install .
# If you want the install to be editable, eg. for development purposes:
# pip install -e .
```

## Usage
```python
from screwrobotics.robot import Robot, Frame
import numpy as np

# Define base poses matrices of all the joints
M0 = np.array([[1, 0,  0, 0],
                [0, 1,  0, 0],
                [0, 0, 1, 0],
                [0, 0,  0, 1]])
M1 = np.array([[1, 0,  0, 0.5],
                [0, 1,  0, 0],
                [0, 0, 1, 0],
                [0, 0,  0, 1]])
M2 = np.array([[1, 0,  0, 1.5],
                [0, 1,  0, 0],
                [0, 0, 1, 0],
                [0, 0,  0, 1]])
M3 = np.array([[1, 0,  0, 2.5],
                [0, 1,  0, 0],
                [0, 0, 1, 0],
                [0, 0,  0, 1]])
M4 = np.array([[1, 0,  0, 3.0],
                [0, 1,  0, 0],
                [0, 0, 1, 0],
                [0, 0,  0, 1]])
M_list = [M0, M1, M2, M3, M4]

# Screw axes for each joint
S1 = np.array([0, 0, 1.0,  0, 0.0, 0])
S2 = np.array([0, 0, 1.0,  0, -1.0, 0])
S3 = np.array([0, 0, 1.0,  0, -2.0, 0])
screw_list = [S1, S2, S3]

# Positions to find forward kinematics for
theta_list = [0, pi / 2, pi / 2]

# Define robot object
r = Robot(M_list, screw_list, Frame.SPACE_FRAME)
# Run forward kinematics and inverse kinematics
fk = r.forward_kinematics([0.1, 0.2, 0.3])
r.inverse_kinematics(fk)
```