# Simulation Results & Visualizations

## 1. Galaxy Positions
Comparison of spatial distributions for different galaxy morphologies.
Yellow cross represents the galactic centre. For elliptical galaxies different colors represent the radial coordinates, while for spiral galaxies they show the populations of different arms 

| Morphology | 1 Million Stars | 100 Million Stars |
| :--- | :---: | :---: |
| **Elliptical** | <img width="1772" height="915" alt="gal_code_1M_elliptical" src="https://github.com/user-attachments/assets/afc3d5d1-221b-49bd-a16a-cd3287ffa777" /> | <img width="1772" height="915" alt="gal_code_100M_elliptical" src="https://github.com/user-attachments/assets/bf43c2aa-c63f-436c-b20c-275dfece5382" /> |
| **Spiral** | <img width="1772" height="915" alt="gal_code_1M" src="https://github.com/user-attachments/assets/5aac0561-2e24-4b1d-b743-5af60340dd01" /> | <img width="1772" height="915" alt="gal_code_100M" src="https://github.com/user-attachments/assets/d4bb1fd3-9f55-44b4-a206-2d4a7377320d" /> |

> **Observation:** Note that the scatter points shrink in size as the stellar sample grows for visualization purposes

## 2. Hertzsprung-Russell (HR) Diagrams
Comparison of stellar evolution tracks using SSE (Single Star Evolution), BSE (Binary Star Evolution) and MESA codes. 

The mass bins considered are [0.11, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100] $M_\odot$
The logP bins considered are [0.15, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5]
The q bins considered are [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

**SSE vs BSE: not that the BSE plots contain more data points and are more heterogeneous because they contain binaries, i.e. twice the number of stars, and because the q ratio offers more stellar tracks on the HR diagram, on top of all the binary products**

| Evolution Code | 1 Million Stars | 100 Million Stars |
| :--- | :---: | :---: |
| **SSE** | <img width="1961" height="978" alt="HR_SSE_1M" src="https://github.com/user-attachments/assets/aff3a5cb-d2ce-4c4c-93c0-1114ae40ead2" /> | <img width="1961" height="978" alt="HR_SSE_100M" src="https://github.com/user-attachments/assets/255a4160-8078-4563-80b6-5caf4b1c6622" />|
| **BSE** | <img width="1961" height="978" alt="HR_BSE_1M" src="https://github.com/user-attachments/assets/ce0dd5aa-cf74-42ce-af57-b7d62e8e2008" /> | <img width="1961" height="978" alt="HR_BSE_100M" src="https://github.com/user-attachments/assets/5a267af3-5019-4f8b-9728-d0a986976aa4" /> |
| **MESA** | <img width="1961" height="978" alt="HR_MESA_1M" src="https://github.com/user-attachments/assets/99111db7-d2e8-4282-b8df-197c201d2585" /> | <img width="1961" height="978" alt="HR_MESA_1M" src="https://github.com/user-attachments/assets/17ec9a35-4242-4ec3-aaf9-541768c80b2d" /> |
