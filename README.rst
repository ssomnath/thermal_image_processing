Thermal Image Processing
=========================
**Suhas Somnath**

This simple standalone Igor-Po procedure file facilitates simple routines for processing results from AFM topography imaging using silicon self-heating AFM cantilevers. It was mainly used for purposes such as:

* Characterizing sensitivity and resolution for thermal topography imaging performed on standardized square grating samples
* Scaling voltage values to nanometers once the sensitivity and resolution were known
* Inverting images when images were acquiring using resistance control
* Adding an offset to all values in a channel - great for centering at 0 V
* Simple edge filtering
* Flipping the image vertically
* Clipping artifacts (spikes) in the signals that would occur at sharp vertical edges of gratings by coercing the values to lie within a min and max
* Replacing the contents of a layer with data present in a text file (relevant if the data of interest was acquired using an alternate source such as a National Instruments Data Acquisition system but if the user prefers to use Asylum Research's / Igor Pro for further image analysis)
* Deconvolving the shape of the cantilever heater from the measured image. See this paper for more information.

Disclaimer
----------
This code was developed from 2010 to 2012 for Asylum Research software versions of 090909. It may or may not work now

Journal papers using this software
-----------------------------------
1. Somnath, Suhas, Elise A. Corbin, and William P. King. "Improved nanotopography sensing via temperature control of a heated atomic force microscope cantilever." IEEE Sensors Journal 11, no. 11 (2011): 2664-2670.
2. Kim, Hoe Joon, Nicolaie Moldovan, Jonathan R. Felts, Suhas Somnath, Zhenting Dai, Tevis DB Jacobs, Robert W. Carpick, John A. Carlisle, and William P. King. "Ultrananocrystalline diamond tip integrated onto a heated atomic force microscope cantilever." Nanotechnology 23, no. 49 (2012): 495302.
3. Lee, Byeonghee, Suhas Somnath, and William P. King. "Fast nanotopography imaging using a high speed cantilever with integrated heater–thermometer." Nanotechnology 24, no. 13 (2013): 135501.
4. Somnath, Suhas, and William P. King. "Heated atomic force cantilever closed loop temperature control and application to high speed nanotopography imaging." Sensors and Actuators A: Physical 192 (2013): 27-33.
5. Liu, Joseph O., Suhas Somnath, and William P. King. "Heated atomic force microscope cantilever with high resistivity for improved temperature sensitivity." Sensors and Actuators A: Physical 201 (2013): 141-147.
6. Somnath, Suhas, Hoe Joon Kim, Huan Hu, and William P. King. "Parallel nanoimaging and nanolithography using a heated microcantilever array." Nanotechnology 25, no. 1 (2013): 014001.
7. Seong, Myunghoon, Suhas Somnath, Hoe Joon Kim, and William P. King. "Parallel nanoimaging using an array of 30 heated microcantilevers." RSC Advances 4, no. 47 (2014): 24747-24754.
8. Somnath, Suhas, and William P. King. "An investigation of heat transfer between a microcantilever and a substrate for improved thermal topography imaging." Nanotechnology 25, no. 36 (2014): 365501.
9. Somnath, Suhas, Joseph O. Liu, Mete Bakir, Craig B. Prater, and William P. King. "Multifunctional atomic force microscope cantilevers with Lorentz force actuation and self-heating capability." Nanotechnology 25, no. 39 (2014): 395501.
