# OSIRiS
Python script to obtain the recombination intermediate DNA sequences of any DNA sequene based on a list of serine integrases of interest.

## Installing

OSIRiS code requires Python 2.7 and installation of biopython, please download the module [here](http://biopython.org/wiki/Download) and to install: `pip install biopython`.

To use OSIRiS script, you will also need a csv file containing your integrase sites of interest. This script works only with serine integrases performing irreversible DNA excision and inversion. The four different sites such as attB, attP, attL and attR of each integrase of interest should be provided. A template csv file with sites of some standard integrases can be found here.
Additionnaly, your DNA sequence of interest must be in genbank format and should be in the directory where the intermediate sequences will be generated. 

## Getting Started

For execution of the function OSIRiS of this python script, go to the line 292 of [osiris_script.py](https://github.com/sguiz/OSIRiS/blob/master/osiris_script.py) and execute the OSIRIS_script function as following: 
`OSIRIS_script(directory, input_file, name_int, directory_list_sites)`
- directory: the directory where you want the DNA sequence to be generated and where the input DNA sequence is.
- input_file: name of the input DNA sequence of interest to be recombined.
- name_int: list of name of integrases of interest.
- directory_list_sites: path of the csv file containing the list of integrase sites.

### Example for Boolean logic

To obtain the intermediate recombinase states of a DNA sequence named sequence_of_interest.gb with the integrases: Bxb1, Tp901 and Int5,
please uncomment the lines 287 to 292 of the code and define the inputs of the functions such as:
`name_int=['Bxb1', 'Tp901', 'Int5']\n
directory='/Users/sarahguiziou/Desktop/osiris'\n
input_file='sequence_of_interest'\n
directory_list_sites=directory+'/List_Integrases.csv'`
And execute the function as follow:
`OSIRIS_script(directory, input_file, name_int, directory_list_sites)`

