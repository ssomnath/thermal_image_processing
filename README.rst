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

This code was developed from 2010 to 2012 for Asylum Research software versions of 090909. It may or may not work now
