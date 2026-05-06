Curent:
    - identified jump after 3.5kpc, regardless of precision, after ||V||>0.05 or so it jumps to another "minima"
    - will try to limit movement distance, or refine somehow
        - start by taking smaller steps and maybe checking values in between, to find the best one w newton?

Done:
	- just implemented models/pot.py for all potentials, as well as effective potential
	- puntequil works correctly (checked against matlabcode)
        - func2 and func2jac must also be correct
    	- der1 and der2 seem to be correct
	- just implemented a first version of eq point continuation. seems to work
    - isopotential and isodensity plots successfully implemented, now using a colormap
    instead of a set of level curves. 
	- but there's a "jump" for the middle L3 eq point (or L1, whatever)
	


