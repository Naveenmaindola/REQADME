# README
README!!!

Code file name :: 1st_order_Kuramoto_multiplexing_delay.jl

Before executing the file--
    1. Please make sure you have Julia installed.
    2. If not already installed, please install 'DifferentialEquations', 'DelimitedFiles', 'Random', 'Printf', 'StatsBase' packages and update Julia.


Input :: Before executing the file please set all the network parameters like number of nodes 'N', 'nearest_neighbour', 'inter', 'intra_1' and 'intra_2' layer coupling and the number of solitary states.


To execute --
    1. open terminal and type 'julia filename.jl'.


Output :: Depending on the values you set for transient time 'tr', 'tf' and 'dt', if it works fine, one will get a file naming 'dn=solit_%dn_hetundel_pht.dat' containing the array of [(tf-tr)/dt] * N.


=============================================================================================================================
=============================================================================================================================


Description of figures ::
"""""""""""""""""""""""

All the data files contain just the phases for all time. After leaving the transients one can easily calculate the corresponding frequency.

    Figure 1 --
    `````````
	Contains the schematic diagram.

    Figure 2 --
    `````````
        Shows the phase and frequency diagram for 1-solitary state when 'intra_1=0.5', 'intra_2=3.0', 'inter=1.0'. (Data file :: 1_Solitary.dat).

    Figure 3 --
    `````````
	Contains the examples of 2, 5 and 10-solitary state with the same parameters as figure 2. (Data file :: 2_Solitary.dat, 5_Solitary.dat, 10_Solitary.dat)

    Figure 4 --
    `````````
	By changing the position of delay, we make the incoherent part of the network spacially localized (Switching between Solitary and Chimera state). (Data file :: 20_Solitary.dat, 		   Chimera(20)nodes).dat) 

    Figure 5 --
    `````````
	For each combination of parameters on X and Y axis, every block of color in Figure 5 shows the difference in the frequency of the solitary node than that of the synchronized chunck.

    Figure 6 --
    `````````
	Figure 6(a) is the phase difference between solitary and synchronized nodes and can be calculated from any of the data file. Similarly, from the difference of those nodes, one can calculate \dot{\eta} for figure (b).
