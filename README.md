# MECA0032 - Flow in turbomachines

## Preliminary design exercise : Redesign of the J79 engine compressor

### Academic Year 2024 – 2025

#### Author: Victor Renkin s2306326

## Table of Contents
- [Introduction](#introduction)
- [Requirements](#requirements)
- [Usage](#usage)
- [Project Structure](#project-structure)
- [AI](#ai)


## Introduction


This project focuses on reverse engineering the renowned General Electric J79 engine. Developed in the 1950s, this highly successful afterburning turbojet engine was initially designed to power the Convair B58 bomber, the first supersonic aircraft capable of reaching Mach 2. Later, it powered the iconic Lockheed F104 Starfighter until the late 1980s and, until recently, the F4 Phantom, one of the most successful and enduring multirole fighters in aviation history.

The J79 engine has been featured in various prototypes, proposed briefly as an alternative engine for the F16, and adapted for aero-derivative gas turbines used in power generation. Table 1 below summarizes the cycle parameters and performance for one of the engine's variants at the design point and sea level static conditions:

| Parameter                         | Value                    |
|-----------------------------------|--------------------------|
| Overall Pressure Ratio (OPR)      | **13.5**                    |
| Turbine Inlet Temperature (TIT)   | 987 °C                  |
| Engine Airflow ($\dot{m}$)        | 77 kg/s                 |
| Thrust ($T$)                      | 52.8 kN / 80 kN         |
| Thrust Specific Fuel Consumption (SFC) | 23.8 g/s/kN / 57.5 g/s/kN |

The engine features a single spool configuration, consisting of a **17-stage axial compressor** with variable inlet guide vanes (VIGV) and six variable stator vanes (VSV), powered by a three-stage axial turbine. Its afterburner employs effusion cooling via turbine exhaust gases through the liner, and the nozzle is a purely converging variable outlet. 

The outer diameter of the engine envelope, including the casing and piping, is approximately 90 cm. The compressor's tip radius remains constant, while the hub diameter increases to compensate for rising air density. For this exercise, the compressor **tip diameter** is assumed to be **80 cm**, with a uniform flow path contraction maintaining a consistent mean radius.

## Reqirements
This project requires the packages listed in the `requirements.txt` file. To install these dependencies, ensure you have [Python](https://www.python.org/) and [pip](https://pip.pypa.io/en/stable/) installed, then run the following command:

```bash
pip install -r requirements.txt
```
## Usage

The project is made in 3 part : the design, the 3D design and also the off-design. 

#### Design
The design part can be computated in the main like this

```bash
python src/main.py
```

#### 3D Design

The 3D design part can be computated in the main like this

```bash
src/3D/main.py
```

#### Off-design

The off-design is also on the main to have the area find by the design, but can be taking outside and juste selected the area for time compilation. 

```bash
python src/main.py
```


## Project Structure

### **`figures/`**
Contains the figures of the project, organized as follows:

- **`design/`**: Includes figures related to the design.
- **`ISR/`**: Contains figures related to the 3D design.
- **`off_design/`**: Contains figures related to the off-design:
  - **`design_n/`**: Figures related to the off-design with $n = n^*$.
  - **`design_n_0_9/`**: Figures related to the off-design with $n = 0.9 \cdot n^*$.
  - **`design_n_1_0_5/`**: Figures related to the off-design with $n = 1.05 \cdot n^*$.

### **`src/`**
Contains the source code of the project, organized as follows:

- **`get_class.py`**: Code defining the classes used in the project:
  - **`atmosphere`**: Represents atmospheric conditions at sea level.
  - **`compressor`**: Contains compressor-related data.
  - **`CompressorTriangle`**: Represents the compressor velocity triangle for pre-design when $v_m = \text{constant}$.
  - **`blade_cascade`**: Contains data for cascades used in off-design analysis.
  - **`station_compressor`**: Represents data for each compressor station.

- **`design/`**: Code for the design part:
  - **`diffusion_interpolation.py`**: Implements the interpolation of diffusion to determine thickness.
  - **`pitchline_design.py`**: Contains the pitchline method for compressor design.
  - **`iterative_design.py`**: Iterative design using the pitchline method to explore possible conditions.

- **`3D/`**: Code for the 3D design part:
  - **`isr.py`**: Implements the ISR method.

- **`off_design/`**: Code for the off-design part:
  - **`off_design.py`**: Implements the off-design method.
  - **`interpolation_off_design.py`**: Contains the interpolation code for the off-design method. This includes two interpolations:
    - Interpolation of the diffusion coefficient $D_{eq}$ to determine thickness.
    - Interpolation of $\zeta$ to calculate the derivative $\frac{d\varepsilon}{d\text{aoa}}$.

- **`main.py`**: Main script to run the project.
- **`VizData.py`**: Code for generating visualizations and figures.


## AI
The AI is used to occasionally correct the code and to reformulate sentences from the report.