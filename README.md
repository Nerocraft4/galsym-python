Curent:
    - will try to limit movement distance, or refine somehow
        - start by taking smaller steps and maybe checking values in between, to find the best one w newton?
    - will for sure try to implement Keller Pseudo-Arc Length KPAL algorithm
    - following main.m, still need to do integration in time, hope it works out well?
Done:
    0.6 Completed aproxlineal.py function, seems to work well
    0.5 Refactored most of the main.py file and created many secondary functions
        - val2func's plot is broken tho
    0.4 isopotential and isodensity plots successfully implemented, now using a colormap instead of a set of level curves. 
    0.3 just implemented a first version of eq point continuation. seems to work
	    - but there's a "jump" for the middle L3 eq point (or L1, whatever)
        - identified jump after 3.5kpc, regardless of precision, after ||V||>0.05 or so it jumps to another "minima"
    0.2 puntequil works correctly (checked against matlabcode)
        - func2 and func2jac must also be correct
    	- der1 and der2 seem to be correct
	0.1 just implemented models/pot.py for all potentials, as well as effective potential
	
	


