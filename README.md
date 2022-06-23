# RRM-RNA-dock
## About The Project

This tool is developed to perform anchor-driven fragment-based docking for a single-stranded RNA that is bound to the amino acid(s) of the RNA recognition motif (RRM) throught the stacking interactions.
Current version 0.1 is limited to docking one ssRNA trinucleotide (further: fragment), two nucleotides of which are bound to the different amino acids of RRM.


## Getting Started

### Prerequisites

- Current version of the tool is developed under Linux (it will not work on Windows or Mac unless a virtual Linux machine is installed) 
- Conda should be installed.
- In order to use RRM-RNA-pip, one needs to know a little bit about interacting with a command line.

### Installation

In order to use RRM-RNA-dock one needs to :

#### 1. Install RRM-RNA-dock and set up its environment

1. Download or clone scipts for the pipeline:
    `git clone https://github.com/AnnaKravchenko/RRM-RNA-dock.git`
2. Create conda environment:
    `conda create --name rrdock --file <your-path>/RRM-RNA-dock/pip-requitements.txt`
3. Define the RRDOCK environment variable for the current terminal session with `export RRDOCK=<your path>/RRM-RNA-dock/` 
    or define it permanently by adding that line in your /home/.bashrc
4. Edit ~/RRM-RNA-dock/mdir/config.ini:
    - replase line "mdir = ~/RRM-RNA-dock/git/mdir/" with "mdir = <your path>/RRM-RNA-dock/git/mdir/" 
    - if needed, change the number of CPUs to use during the docking by replasing line "docking_cpu = 8" with "docking_cpu = <your number>"

#### 2. Install ATTRACT docking engine and set up its environment

1. Install ATTRACT from here: https://github.com/sjdv1982/attract 
2. Define the ATTRACTDIR environment variable for the current terminal session with `export ATTRACTDIR=<your path>/attract/tools/` 
    or define it permanently by adding that line in your /home/.bashrc
3. Define the ATTRACTTOOLS environment variable for the current terminal session with `export ATTRACTTOOLS=<your path>/attract/bin/` 
    or define it permanently by adding that line in your /home/.bashrc

#### 3. Download fragments library and set up path to it 

1. Download library from here: https://zenodo.org/record/6483823#.YrQttNJByXJ, you are looking for the file trinucl_clust1A_ATTRACT.tar
2. Follow directions in README to decompress library files from motif.npy to motif.pdb
3. Define the LIBRARY environment variable for the current terminal session with `export LIBRARY=<your path>/fraglib/`
    or define it permanently by adding that line in your /home/.bashrc

### Usage

1. Activate required environment by `conda activate rrdock`
2. - To run pipeline for the protiein Q61474_RRM1 and ssRNA UUU, 
    - where first nucleoide is expected to stack to RNP2, 
    - and second nucleotide is to stack to RNP1, 
    - and to save the results to the folder home/test/, run the next command line:
    `python3 $RRDOCK/pip.py -id Q61474 -rrm 1 -seq UUU -ancNucB1 1 -ancNucB3 2 -wdir home/test/`

## Contributing
Contributions are what makes the open source community such an amazing place to learn, inspire, and create. Any contributions you make are greatly appreciated.

If you have a suggestion that would make this better, please fork the repo and create a pull request. 

1. Fork the Project
2. Create your Feature Branch (git checkout -b feature/YourFeature)
3. Commit your Changes (git commit -m 'Add some YourFeature')
3. Push to the Branch (git push origin feature/YourFeature)
4. Open a Pull Request


## License

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Contacts
Anna Kravchenko - anna.kravchenko@loria.fr 

Hrishikesh Dhondge - hrishikesh.dhondge@loria.fr

Project Link: http://rnact.eu/

## Acknowledgments
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No. 813239.

