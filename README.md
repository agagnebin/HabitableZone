# Habitable Zone
The goal of this project is to create a plot of the orbit of an Earth-like planet in the hospitable zone of any main sequence star as a function of time.

To run this project:

First, compile FinalProject.cpp and then run the executable along with an argument for the spectral class. The accepted arguments are O, B, A, F, G, K, and M (can be lowercase or uppercase). If you enter anything else, it will say that it isn't an accepted spectral class.
This will create a file, "Final_Positions.txt" which will be read and plotted in Python.

To run the Python file:
 python3 FinalProject.py 

To animate the plots:  
cd Plots
ffmpeg -fflags +genpts -r 50 -i fig%04d.png -vcodec mpeg4 -qmax 1000 orbit.mp4


##References
https://arxiv.org/pdf/1301.6674.pdf
https://en.wikipedia.org/wiki/Orbital_speed#Instantaneous_orbital_speed
https://en.wikipedia.org/wiki/Initial_mass_function#Kroupa_(2001)
Classical Mechanics by John R. Taylor
